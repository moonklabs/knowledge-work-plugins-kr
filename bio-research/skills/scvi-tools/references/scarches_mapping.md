# scArches를 이용한 레퍼런스 매핑

본 레퍼런스는 scArches를 사용하여 처음부터 재학습하지 않고 사전 학습된 레퍼런스 모델에 쿼리 데이터를 매핑하는 방법을 다룹니다.

## 개요

scArches(single-cell architecture surgery)는 다음을 가능하게 합니다:
- 기존 레퍼런스 아틀라스에 새 데이터 매핑
- 새로운 배치/연구로 모델 확장
- 전체 재학습 없이 전이 학습 수행
- 쿼리 통합 시 레퍼런스 구조 보존

## scArches 사용 시점

| 시나리오 | 접근법 |
|----------|--------|
| 기존 아틀라스에 쿼리 매핑 | scArches 쿼리 매핑 |
| 새 데이터로 아틀라스 확장 | scArches 모델 수술 |
| 사전 학습된 모델 없음 | scANVI를 처음부터 학습 |
| 쿼리가 레퍼런스와 매우 다름 | 재학습 고려 |

## 사전 요구 사항

```python
import scvi
import scanpy as sc
import numpy as np

print(f"scvi-tools version: {scvi.__version__}")
```

## 워크플로우 1: 사전 학습된 레퍼런스에 쿼리 매핑

### 1단계: 사전 학습된 레퍼런스 모델 로드

```python
# Load saved reference model
# The model must have been trained with scvi-tools
reference_model = scvi.model.SCVI.load("reference_model/")

# Or load scANVI for label transfer
reference_model = scvi.model.SCANVI.load("reference_scanvi_model/")

# Check model info
print(f"Model type: {type(reference_model)}")
print(f"Training data shape: {reference_model.adata.shape}")
```

### 2단계: 쿼리 데이터 준비

```python
# Load query data
adata_query = sc.read_h5ad("query_data.h5ad")

# CRITICAL: Match genes to reference
reference_genes = reference_model.adata.var_names
query_genes = adata_query.var_names

# Check overlap
common_genes = reference_genes.intersection(query_genes)
print(f"Reference genes: {len(reference_genes)}")
print(f"Query genes: {len(query_genes)}")
print(f"Overlap: {len(common_genes)}")

# Subset query to reference genes
adata_query = adata_query[:, reference_genes].copy()

# Handle missing genes (filled with zeros automatically by prepare_query_anndata)
```

### 3단계: 쿼리 AnnData 준비

```python
# Store raw counts
adata_query.layers["counts"] = adata_query.X.copy()

# Prepare query for mapping
# This aligns the query data structure to match the reference
scvi.model.SCVI.prepare_query_anndata(adata_query, reference_model)
```

### 4단계: 쿼리 모델 생성

```python
# Create query model from reference
# This initializes with reference weights
query_model = scvi.model.SCVI.load_query_data(
    adata_query,
    reference_model
)

# The query model inherits:
# - Reference architecture
# - Reference encoder weights (frozen by default)
# - Decoder is fine-tuned for query
```

### 5단계: 쿼리에 대한 미세 조정

```python
# Fine-tune the query model
# This adjusts decoder weights for query-specific effects
query_model.train(
    max_epochs=200,
    plan_kwargs={
        "weight_decay": 0.0  # Less regularization for fine-tuning
    }
)

# Check training
query_model.history['elbo_train'].plot()
```

### 6단계: 쿼리 표현 획득

```python
# Get latent representation
# Query cells are embedded in same space as reference
adata_query.obsm["X_scVI"] = query_model.get_latent_representation()

# Visualize
sc.pp.neighbors(adata_query, use_rep="X_scVI")
sc.tl.umap(adata_query)
sc.pl.umap(adata_query, color=['cell_type', 'batch'])
```

## 워크플로우 2: 레이블 전이를 포함한 scANVI 쿼리 매핑

레퍼런스에서 쿼리로 세포 유형 레이블을 전이하는 방법입니다:

### 1단계: scANVI 레퍼런스 로드

```python
# Reference must be scANVI model (trained with labels)
reference_scanvi = scvi.model.SCANVI.load("scanvi_reference/")

# Check available labels
print("Reference cell types:")
print(reference_scanvi.adata.obs['cell_type'].value_counts())
```

### 2단계: 쿼리 준비 및 매핑

```python
# Prepare query
adata_query.layers["counts"] = adata_query.X.copy()
adata_query = adata_query[:, reference_scanvi.adata.var_names].copy()

scvi.model.SCANVI.prepare_query_anndata(adata_query, reference_scanvi)

# Create query model
query_scanvi = scvi.model.SCANVI.load_query_data(
    adata_query,
    reference_scanvi
)

# Fine-tune
query_scanvi.train(
    max_epochs=100,
    plan_kwargs={"weight_decay": 0.0}
)
```

### 3단계: 예측 결과 획득

```python
# Predict cell types
predictions = query_scanvi.predict()
adata_query.obs["predicted_cell_type"] = predictions

# Get prediction probabilities
soft_predictions = query_scanvi.predict(soft=True)
adata_query.obs["prediction_confidence"] = soft_predictions.max(axis=1)

# Latent representation
adata_query.obsm["X_scANVI"] = query_scanvi.get_latent_representation()

# Visualize predictions
sc.pp.neighbors(adata_query, use_rep="X_scANVI")
sc.tl.umap(adata_query)
sc.pl.umap(adata_query, color=['predicted_cell_type', 'prediction_confidence'])
```

### 4단계: 예측 평가

```python
# Distribution of predictions
print(adata_query.obs['predicted_cell_type'].value_counts())

# Confidence statistics
print(f"Mean confidence: {adata_query.obs['prediction_confidence'].mean():.3f}")
print(f"Low confidence (<0.5): {(adata_query.obs['prediction_confidence'] < 0.5).sum()}")

# Filter low-confidence predictions
high_conf = adata_query[adata_query.obs['prediction_confidence'] >= 0.7].copy()
print(f"High confidence cells: {len(high_conf)} ({len(high_conf)/len(adata_query)*100:.1f}%)")
```

## 워크플로우 3: 모델 수술 (레퍼런스 확장)

기존 레퍼런스 모델을 새 데이터로 확장합니다:

### 1단계: 레퍼런스 레이어 고정

```python
# Load reference model
reference_model = scvi.model.SCVI.load("reference_model/")

# Get reference representation (before surgery)
adata_ref = reference_model.adata
adata_ref.obsm["X_scVI_before"] = reference_model.get_latent_representation()
```

### 2단계: 결합 데이터 준비

```python
# Add batch information
adata_ref.obs["dataset"] = "reference"
adata_query.obs["dataset"] = "query"

# Combine
adata_combined = sc.concat([adata_ref, adata_query])
adata_combined.layers["counts"] = adata_combined.X.copy()
```

### 3단계: 수술 접근법

```python
# Option A: Use load_query_data (recommended)
scvi.model.SCVI.prepare_query_anndata(adata_query, reference_model)
extended_model = scvi.model.SCVI.load_query_data(adata_query, reference_model)
extended_model.train(max_epochs=200)

# Option B: Retrain with combined data (if query is large)
# This doesn't preserve reference exactly but may give better results
scvi.model.SCVI.setup_anndata(
    adata_combined,
    layer="counts",
    batch_key="dataset"
)
new_model = scvi.model.SCVI(adata_combined, n_latent=30)
new_model.train(max_epochs=200)
```

## 결합 시각화

레퍼런스와 쿼리를 함께 시각화합니다:

```python
# Get latent representations
adata_ref.obsm["X_scVI"] = reference_model.get_latent_representation()
adata_query.obsm["X_scVI"] = query_model.get_latent_representation()

# Combine for visualization
adata_ref.obs["source"] = "reference"
adata_query.obs["source"] = "query"
adata_combined = sc.concat([adata_ref, adata_query])

# Compute joint UMAP
sc.pp.neighbors(adata_combined, use_rep="X_scVI")
sc.tl.umap(adata_combined)

# Visualize
import matplotlib.pyplot as plt
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

sc.pl.umap(adata_combined, color="source", ax=axes[0], show=False, title="Source")
sc.pl.umap(adata_combined, color="cell_type", ax=axes[1], show=False, title="Cell Type")
sc.pl.umap(adata_combined, color="batch", ax=axes[2], show=False, title="Batch")

plt.tight_layout()
```

## 공개 아틀라스 모델 사용

### HuggingFace Model Hub에서 가져오기

```python
from huggingface_hub import hf_hub_download

# Download model files
model_dir = hf_hub_download(
    repo_id="scvi-tools/model-name",  # Replace with actual repo
    filename="model.pt",
    local_dir="./downloaded_model/"
)

# Load model
atlas_model = scvi.model.SCANVI.load(model_dir)
```

### CellxGene에서 가져오기

```python
# Many CellxGene datasets provide pre-trained models
# Check dataset documentation for model availability
# https://cellxgene.cziscience.com/

# Example workflow:
# 1. Download reference dataset and model
# 2. Load model: model = scvi.model.SCANVI.load("cellxgene_model/")
# 3. Map your query data using steps above
```

## 전체 파이프라인

```python
def map_query_to_reference(
    adata_query,
    reference_model_path,
    model_type="scanvi",
    max_epochs=100,
    confidence_threshold=0.5
):
    """
    Map query data to pre-trained reference model.

    Parameters
    ----------
    adata_query : AnnData
        Query data with raw counts
    reference_model_path : str
        Path to saved reference model
    model_type : str
        "scvi" or "scanvi"
    max_epochs : int
        Fine-tuning epochs
    confidence_threshold : float
        Minimum prediction confidence (for scANVI)

    Returns
    -------
    Mapped AnnData with predictions (if scANVI)
    """
    import scvi

    # Load reference
    if model_type == "scanvi":
        reference_model = scvi.model.SCANVI.load(reference_model_path)
        ModelClass = scvi.model.SCANVI
    else:
        reference_model = scvi.model.SCVI.load(reference_model_path)
        ModelClass = scvi.model.SCVI

    # Prepare query
    adata_query = adata_query.copy()
    adata_query = adata_query[:, reference_model.adata.var_names].copy()
    adata_query.layers["counts"] = adata_query.X.copy()

    # Map query
    ModelClass.prepare_query_anndata(adata_query, reference_model)
    query_model = ModelClass.load_query_data(adata_query, reference_model)

    # Fine-tune
    query_model.train(
        max_epochs=max_epochs,
        plan_kwargs={"weight_decay": 0.0}
    )

    # Get results
    rep_key = "X_scANVI" if model_type == "scanvi" else "X_scVI"
    adata_query.obsm[rep_key] = query_model.get_latent_representation()

    if model_type == "scanvi":
        adata_query.obs["predicted_cell_type"] = query_model.predict()
        soft = query_model.predict(soft=True)
        adata_query.obs["prediction_confidence"] = soft.max(axis=1)
        adata_query.obs["confident"] = adata_query.obs["prediction_confidence"] >= confidence_threshold

    # Compute UMAP
    sc.pp.neighbors(adata_query, use_rep=rep_key)
    sc.tl.umap(adata_query)

    return adata_query, query_model


# Usage
adata_mapped, model = map_query_to_reference(
    adata_query,
    "reference_scanvi_model/",
    model_type="scanvi"
)

# Visualize
sc.pl.umap(adata_mapped, color=['predicted_cell_type', 'prediction_confidence'])
```

## 문제 해결

| 문제 | 원인 | 해결 방법 |
|------|------|-----------|
| 유전자 불일치 | 유전자 명명법 차이 | 유전자 ID 변환 (Ensembl ↔ Symbol) |
| 낮은 신뢰도 다수 | 쿼리에 새로운 세포 유형 존재 | 낮은 신뢰도 세포를 수동 주석 처리 |
| 매핑 품질 저하 | 쿼리와 레퍼런스 차이가 큼 | 결합 데이터로 재학습 고려 |
| 메모리 오류 | 대규모 쿼리 | 배치 단위로 처리 |
| 버전 불일치 | scvi-tools 버전 차이 | 레퍼런스 학습 시 사용된 동일 버전 사용 |

## 주요 참고문헌

- Lotfollahi et al. (2022) "Mapping single-cell data to reference atlases by transfer learning"
- Xu et al. (2021) "Probabilistic harmonization and annotation of single-cell transcriptomics data with deep generative models"
