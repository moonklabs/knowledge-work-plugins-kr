# scANVI를 이용한 라벨 전이 및 레퍼런스 매핑

이 레퍼런스는 scANVI를 사용하여 레퍼런스 아틀라스에서 쿼리 데이터로 세포 유형 주석을 전이하는 방법을 다룹니다.

## 개요

레퍼런스 매핑("라벨 전이"라고도 함)은 주석이 달린 레퍼런스 데이터에서 사전 학습된 모델을 사용하여 새로운 미주석 쿼리 데이터의 세포 유형을 예측합니다. 이 방법은 재클러스터링보다 빠르고 연구 간 더 일관된 결과를 제공합니다.

scANVI가 이 작업에 뛰어난 이유:
- 레퍼런스와 쿼리를 공유 공간에 동시 임베딩
- 확률적 라벨 전이
- 레퍼런스와 쿼리 간 배치 효과 처리

## 레퍼런스 매핑을 사용해야 하는 경우

- 기존 아틀라스를 사용하여 새 데이터셋에 주석 부여
- 여러 연구에 걸쳐 일관된 주석 수행
- 속도: 재클러스터링 및 수동 주석 불필요
- 품질: 전문가가 큐레이션한 레퍼런스 주석 활용

## 워크플로우 옵션

1. **새 모델 학습**: 레퍼런스에서 scANVI를 학습한 후 쿼리 매핑
2. **사전 학습 모델 사용**: 기존 모델 로드 (예: Model Hub에서)
3. **scArches**: 쿼리 데이터로 기존 모델 확장 (레퍼런스 보존)

## 옵션 1: 레퍼런스에서 scANVI 학습

### 1단계: 레퍼런스 데이터 준비

```python
import scvi
import scanpy as sc

# Load reference atlas
adata_ref = sc.read_h5ad("reference_atlas.h5ad")

# Check annotations
print(f"Reference cells: {adata_ref.n_obs}")
print(f"Cell types: {adata_ref.obs['cell_type'].nunique()}")
print(adata_ref.obs['cell_type'].value_counts())

# Ensure raw counts
adata_ref.layers["counts"] = adata_ref.raw.X.copy() if adata_ref.raw else adata_ref.X.copy()

# HVG selection
sc.pp.highly_variable_genes(
    adata_ref,
    n_top_genes=3000,
    flavor="seurat_v3",
    batch_key="batch" if "batch" in adata_ref.obs else None,
    layer="counts"
)
adata_ref = adata_ref[:, adata_ref.var["highly_variable"]].copy()
```

### 2단계: 레퍼런스에서 scANVI 학습

```python
# First train scVI (unlabeled)
scvi.model.SCVI.setup_anndata(
    adata_ref,
    layer="counts",
    batch_key="batch"
)

scvi_ref = scvi.model.SCVI(adata_ref, n_latent=30)
scvi_ref.train(max_epochs=200)

# Initialize scANVI from scVI
scanvi_ref = scvi.model.SCANVI.from_scvi_model(
    scvi_ref,
    labels_key="cell_type",
    unlabeled_category="Unknown"
)

# Train scANVI
scanvi_ref.train(max_epochs=50)

# Save for later use
scanvi_ref.save("scanvi_reference_model/")
```

### 3단계: 쿼리 데이터 준비

```python
# Load query data
adata_query = sc.read_h5ad("query_data.h5ad")

# CRITICAL: Use same genes as reference
common_genes = adata_ref.var_names.intersection(adata_query.var_names)
print(f"Common genes: {len(common_genes)}")

# Subset query to reference genes
adata_query = adata_query[:, adata_ref.var_names].copy()

# Handle missing genes (set to 0)
missing_genes = set(adata_ref.var_names) - set(adata_query.var_names)
if missing_genes:
    # Add missing genes with zero expression
    import numpy as np
    from scipy.sparse import csr_matrix

    zero_matrix = csr_matrix((adata_query.n_obs, len(missing_genes)))
    # ... concat and reorder to match reference

# Store counts
adata_query.layers["counts"] = adata_query.X.copy()
```

### 4단계: 쿼리를 레퍼런스에 매핑

```python
# Prepare query data for mapping
scvi.model.SCANVI.prepare_query_anndata(adata_query, scanvi_ref)

# Create query model from reference
scanvi_query = scvi.model.SCANVI.load_query_data(
    adata_query,
    scanvi_ref
)

# Fine-tune on query (optional but recommended)
scanvi_query.train(
    max_epochs=100,
    plan_kwargs={"weight_decay": 0.0}
)

# Get predictions
adata_query.obs["predicted_cell_type"] = scanvi_query.predict()

# Get prediction probabilities
soft_predictions = scanvi_query.predict(soft=True)
adata_query.obs["prediction_score"] = soft_predictions.max(axis=1)
```

### 5단계: 예측 결과 평가

```python
# Confidence scores
print(f"Mean prediction confidence: {adata_query.obs['prediction_score'].mean():.3f}")

# Low confidence predictions
low_conf = adata_query.obs['prediction_score'] < 0.5
print(f"Low confidence cells: {low_conf.sum()} ({low_conf.mean()*100:.1f}%)")

# Visualize
sc.pp.neighbors(adata_query, use_rep="X_scANVI")
sc.tl.umap(adata_query)
sc.pl.umap(adata_query, color=['predicted_cell_type', 'prediction_score'])
```

## 옵션 2: 사전 학습 모델 사용

### Model Hub에서 로드

```python
# scvi-tools maintains models on HuggingFace
# Check: https://huggingface.co/scvi-tools

# Example: Load pre-trained model
from huggingface_hub import hf_hub_download

model_path = hf_hub_download(
    repo_id="scvi-tools/example-model",
    filename="model.pt"
)

# Load model
model = scvi.model.SCANVI.load(model_path, adata=adata_query)
```

### 출판된 아틀라스에서 로드

```python
# Many atlases provide pre-trained models
# Example workflow with CellTypist-style model

# Download reference model
# model = scvi.model.SCANVI.load("atlas_model/", adata=adata_query)
```

## 옵션 3: 점진적 업데이트를 위한 scArches

scArches는 처음부터 재학습하지 않고 레퍼런스 모델을 확장합니다:

```python
# Load existing reference model
scanvi_ref = scvi.model.SCANVI.load("reference_model/")

# Surgery: prepare for query integration
scanvi_ref.freeze_layers()

# Map query data
scvi.model.SCANVI.prepare_query_anndata(adata_query, scanvi_ref)
scanvi_query = scvi.model.SCANVI.load_query_data(adata_query, scanvi_ref)

# Train only query-specific parameters
scanvi_query.train(
    max_epochs=200,
    plan_kwargs={"weight_decay": 0.0}
)
```

## 레퍼런스와 쿼리의 동시 시각화

```python
# Concatenate for joint visualization
adata_ref.obs["dataset"] = "reference"
adata_query.obs["dataset"] = "query"

# Get latent representations
adata_ref.obsm["X_scANVI"] = scanvi_ref.get_latent_representation()
adata_query.obsm["X_scANVI"] = scanvi_query.get_latent_representation()

# Combine
adata_combined = sc.concat([adata_ref, adata_query])

# Compute combined UMAP
sc.pp.neighbors(adata_combined, use_rep="X_scANVI")
sc.tl.umap(adata_combined)

# Plot
sc.pl.umap(
    adata_combined,
    color=["dataset", "cell_type", "predicted_cell_type"],
    ncols=2
)
```

## 예측 품질 관리

### 신뢰도 필터링

```python
# Filter predictions by confidence
confidence_threshold = 0.7

high_conf = adata_query[adata_query.obs['prediction_score'] >= confidence_threshold].copy()
low_conf = adata_query[adata_query.obs['prediction_score'] < confidence_threshold].copy()

print(f"High confidence: {len(high_conf)} ({len(high_conf)/len(adata_query)*100:.1f}%)")
print(f"Low confidence: {len(low_conf)} ({len(low_conf)/len(adata_query)*100:.1f}%)")
```

### 마커 유전자 검증

```python
# Validate predictions with known markers
markers = {
    'T cells': ['CD3D', 'CD3E'],
    'B cells': ['CD19', 'MS4A1'],
    'Monocytes': ['CD14', 'LYZ']
}

for ct, genes in markers.items():
    ct_cells = adata_query[adata_query.obs['predicted_cell_type'] == ct]
    if len(ct_cells) > 0:
        for gene in genes:
            if gene in adata_query.var_names:
                expr = ct_cells[:, gene].X.mean()
                print(f"{ct} - {gene}: {expr:.3f}")
```

## 전체 파이프라인

```python
def transfer_labels(
    adata_ref,
    adata_query,
    cell_type_key="cell_type",
    batch_key=None,
    n_top_genes=3000,
    confidence_threshold=0.5
):
    """
    Transfer cell type labels from reference to query.

    Parameters
    ----------
    adata_ref : AnnData
        Annotated reference data
    adata_query : AnnData
        Unannotated query data
    cell_type_key : str
        Column with cell type annotations in reference
    batch_key : str, optional
        Batch column
    n_top_genes : int
        Number of HVGs
    confidence_threshold : float
        Minimum confidence for predictions

    Returns
    -------
    AnnData with predictions
    """
    import scvi
    import scanpy as sc

    # Prepare reference
    adata_ref = adata_ref.copy()
    adata_ref.layers["counts"] = adata_ref.X.copy()

    sc.pp.highly_variable_genes(
        adata_ref,
        n_top_genes=n_top_genes,
        flavor="seurat_v3",
        batch_key=batch_key,
        layer="counts"
    )
    adata_ref = adata_ref[:, adata_ref.var["highly_variable"]].copy()

    # Train reference model
    scvi.model.SCVI.setup_anndata(adata_ref, layer="counts", batch_key=batch_key)
    scvi_ref = scvi.model.SCVI(adata_ref, n_latent=30)
    scvi_ref.train(max_epochs=200)

    scanvi_ref = scvi.model.SCANVI.from_scvi_model(
        scvi_ref,
        labels_key=cell_type_key,
        unlabeled_category="Unknown"
    )
    scanvi_ref.train(max_epochs=50)

    # Prepare query
    adata_query = adata_query[:, adata_ref.var_names].copy()
    adata_query.layers["counts"] = adata_query.X.copy()

    # Map query
    scvi.model.SCANVI.prepare_query_anndata(adata_query, scanvi_ref)
    scanvi_query = scvi.model.SCANVI.load_query_data(adata_query, scanvi_ref)
    scanvi_query.train(max_epochs=100, plan_kwargs={"weight_decay": 0.0})

    # Get predictions
    adata_query.obs["predicted_cell_type"] = scanvi_query.predict()
    soft = scanvi_query.predict(soft=True)
    adata_query.obs["prediction_score"] = soft.max(axis=1)

    # Mark low confidence
    adata_query.obs["confident_prediction"] = adata_query.obs["prediction_score"] >= confidence_threshold

    # Add latent representation
    adata_query.obsm["X_scANVI"] = scanvi_query.get_latent_representation()

    return adata_query, scanvi_ref, scanvi_query

# Usage
adata_annotated, ref_model, query_model = transfer_labels(
    adata_ref,
    adata_query,
    cell_type_key="cell_type"
)

# Visualize
sc.pp.neighbors(adata_annotated, use_rep="X_scANVI")
sc.tl.umap(adata_annotated)
sc.pl.umap(adata_annotated, color=['predicted_cell_type', 'prediction_score'])
```

## 문제 해결

| 문제 | 원인 | 해결 방법 |
|------|------|----------|
| 낮은 신뢰도 예측이 많음 | 쿼리에 새로운 세포 유형 존재 | 낮은 신뢰도 세포를 수동 주석 |
| 잘못된 예측 | 레퍼런스가 조직과 불일치 | 조직에 적합한 레퍼런스 사용 |
| 유전자 불일치 | 다른 유전자 명명 체계 | 유전자 ID 변환 |
| 모두 동일한 예측 | 쿼리가 너무 다름 | 데이터 품질 확인, 다른 레퍼런스 시도 |

## 주요 참고문헌

- Xu et al. (2021) "Probabilistic harmonization and annotation of single-cell transcriptomics data with deep generative models"
- Lotfollahi et al. (2022) "Mapping single-cell data to reference atlases by transfer learning"
