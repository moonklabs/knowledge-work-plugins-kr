# scVI 및 scANVI를 이용한 단일세포 RNA 시퀀싱 데이터 통합

본 레퍼런스는 scVI(비지도 학습) 및 scANVI(세포 유형 레이블을 활용한 준지도 학습)를 사용한 배치 보정 및 데이터셋 통합을 다룹니다.

## 개요

단일세포 데이터셋에는 다음과 같은 원인으로 배치 효과가 자주 발생합니다:
- 서로 다른 기증자/환자
- 서로 다른 실험 배치
- 서로 다른 기술 (10x v2 vs v3)
- 서로 다른 연구

scVI와 scANVI는 배치 효과를 제거하면서 생물학적 변이를 보존하는 공유 잠재 공간을 학습합니다.

## 모델 선택 기준

| 모델 | 사용 시점 | 레이블 필요 여부 |
|------|----------|-----------------|
| **scVI** | 레이블 없이 탐색적 분석 수행 시 | 아니오 |
| **scANVI** | 부분/전체 레이블이 있고 더 나은 보존 필요 시 | 예 (부분 레이블 가능) |

## scVI 통합 워크플로우

### 1단계: 데이터 준비

```python
import scvi
import scanpy as sc

# Load datasets
adata1 = sc.read_h5ad("dataset1.h5ad")
adata2 = sc.read_h5ad("dataset2.h5ad")

# Add batch annotation
adata1.obs["batch"] = "batch1"
adata2.obs["batch"] = "batch2"

# Concatenate
adata = sc.concat([adata1, adata2], label="batch")

# Ensure we have raw counts
# If data is normalized, recover from .raw
if hasattr(adata, 'raw') and adata.raw is not None:
    adata = adata.raw.to_adata()

# Store counts
adata.layers["counts"] = adata.X.copy()
```

### 2단계: 배치 간 고변이 유전자(HVG) 선택

```python
# Select HVGs considering batch
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    flavor="seurat_v3",
    batch_key="batch",
    layer="counts"
)

# Subset to HVGs
adata = adata[:, adata.var["highly_variable"]].copy()
```

### 3단계: scVI 설정 및 학습

```python
# Register data with scVI
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="batch"
)

# Create model
model = scvi.model.SCVI(
    adata,
    n_latent=30,          # Latent dimensions
    n_layers=2,           # Encoder/decoder depth
    gene_likelihood="nb"  # negative binomial (or "zinb")
)

# Train
model.train(
    max_epochs=200,
    early_stopping=True,
    early_stopping_patience=10,
    batch_size=128
)

# Plot training history
model.history["elbo_train"].plot()
```

### 4단계: 통합된 표현 획득

```python
# Get latent representation
adata.obsm["X_scVI"] = model.get_latent_representation()

# Use for clustering and visualization
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=15)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=1.0)

# Visualize integration
sc.pl.umap(adata, color=["batch", "leiden"], ncols=2)
```

### 5단계: 모델 저장

```python
# Save model for later use
model.save("scvi_model/")

# Load model
model = scvi.model.SCVI.load("scvi_model/", adata=adata)
```

## scANVI 통합 워크플로우

scANVI는 scVI에 세포 유형 레이블을 추가하여 생물학적 보존을 향상시킵니다.

### 1단계: 레이블이 포함된 데이터 준비

```python
# Labels should be in adata.obs
# Use "Unknown" for unlabeled cells
print(adata.obs["cell_type"].value_counts())

# For partially labeled data
# Mark unlabeled cells
adata.obs["cell_type_scanvi"] = adata.obs["cell_type"].copy()
# adata.obs.loc[unlabeled_mask, "cell_type_scanvi"] = "Unknown"
```

### 2단계: 옵션 A - scANVI를 처음부터 학습

```python
# Setup for scANVI
scvi.model.SCANVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="batch",
    labels_key="cell_type"
)

# Create model
scanvi_model = scvi.model.SCANVI(
    adata,
    n_latent=30,
    n_layers=2
)

# Train
scanvi_model.train(max_epochs=200)
```

### 2단계: 옵션 B - scVI에서 scANVI 초기화 (권장)

```python
# First train scVI
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
scvi_model = scvi.model.SCVI(adata, n_latent=30)
scvi_model.train(max_epochs=200)

# Initialize scANVI from scVI
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    labels_key="cell_type",
    unlabeled_category="Unknown"  # For partially labeled data
)

# Fine-tune scANVI (fewer epochs needed)
scanvi_model.train(max_epochs=50)
```

### 3단계: 결과 획득

```python
# Latent representation
adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation()

# Predicted labels for unlabeled cells
predictions = scanvi_model.predict()
adata.obs["predicted_cell_type"] = predictions

# Prediction probabilities
soft_predictions = scanvi_model.predict(soft=True)

# Visualization
sc.pp.neighbors(adata, use_rep="X_scANVI")
sc.tl.umap(adata)
sc.pl.umap(adata, color=["batch", "cell_type", "predicted_cell_type"])
```

## 통합 품질 비교

### 시각적 평가

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Before integration (on PCA)
sc.pp.pca(adata)
sc.pl.pca(adata, color="batch", ax=axes[0], title="Before (PCA)", show=False)

# After scVI
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata)
sc.pl.umap(adata, color="batch", ax=axes[1], title="After scVI", show=False)

# After scANVI
sc.pp.neighbors(adata, use_rep="X_scANVI")
sc.tl.umap(adata)
sc.pl.umap(adata, color="batch", ax=axes[2], title="After scANVI", show=False)

plt.tight_layout()
```

### 정량적 지표 (scib)

```python
# pip install scib-metrics

from scib_metrics.benchmark import Benchmarker

bm = Benchmarker(
    adata,
    batch_key="batch",
    label_key="cell_type",
    embedding_obsm_keys=["X_pca", "X_scVI", "X_scANVI"]
)

bm.benchmark()
bm.plot_results_table()
```

## 차등 발현 분석

scVI는 배치 효과를 고려한 차등 발현 분석을 제공합니다:

```python
# DE between groups
de_results = model.differential_expression(
    groupby="cell_type",
    group1="T cells",
    group2="B cells"
)

# Filter significant
de_sig = de_results[
    (de_results["is_de_fdr_0.05"] == True) &
    (abs(de_results["lfc_mean"]) > 1)
]

print(de_sig.head(20))
```

## 고급: 다중 범주형 공변량

```python
# Include additional covariates beyond batch
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="batch",
    categorical_covariate_keys=["donor", "technology"]
)

model = scvi.model.SCVI(adata, n_latent=30)
model.train()
```

## 학습 팁

### 대규모 데이터셋 (>100k 세포)

```python
model.train(
    max_epochs=100,      # Fewer epochs needed
    batch_size=256,      # Larger batches
    train_size=0.9,      # Less validation
    early_stopping=True
)
```

### 소규모 데이터셋 (<10k 세포)

```python
model = scvi.model.SCVI(
    adata,
    n_latent=10,         # Smaller latent space
    n_layers=1,          # Simpler model
    dropout_rate=0.2     # More regularization
)

model.train(
    max_epochs=400,
    batch_size=64
)
```

### 학습 모니터링

```python
# Check training curves
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
ax.plot(model.history["elbo_train"], label="Train")
ax.plot(model.history["elbo_validation"], label="Validation")
ax.set_xlabel("Epoch")
ax.set_ylabel("ELBO")
ax.legend()

# Should see convergence without overfitting
```

## 전체 파이프라인

```python
def integrate_datasets(
    adatas,
    batch_key="batch",
    labels_key=None,
    n_top_genes=2000,
    n_latent=30
):
    """
    Integrate multiple scRNA-seq datasets.

    Parameters
    ----------
    adatas : dict
        Dictionary of {batch_name: AnnData}
    batch_key : str
        Key for batch annotation
    labels_key : str, optional
        Key for cell type labels (uses scANVI if provided)
    n_top_genes : int
        Number of HVGs
    n_latent : int
        Latent dimensions

    Returns
    -------
    AnnData with integrated representation
    """
    import scvi
    import scanpy as sc

    # Add batch labels and concatenate
    for batch_name, adata in adatas.items():
        adata.obs[batch_key] = batch_name

    adata = sc.concat(list(adatas.values()), label=batch_key)

    # Store counts
    adata.layers["counts"] = adata.X.copy()

    # HVG selection
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor="seurat_v3",
        batch_key=batch_key,
        layer="counts"
    )
    adata = adata[:, adata.var["highly_variable"]].copy()

    # Train model
    if labels_key and labels_key in adata.obs.columns:
        # Use scANVI
        scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=batch_key)
        scvi_model = scvi.model.SCVI(adata, n_latent=n_latent)
        scvi_model.train(max_epochs=200)

        model = scvi.model.SCANVI.from_scvi_model(
            scvi_model,
            labels_key=labels_key,
            unlabeled_category="Unknown"
        )
        model.train(max_epochs=50)
        rep_key = "X_scANVI"
    else:
        # Use scVI
        scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=batch_key)
        model = scvi.model.SCVI(adata, n_latent=n_latent)
        model.train(max_epochs=200)
        rep_key = "X_scVI"

    # Add representation
    adata.obsm[rep_key] = model.get_latent_representation()

    # Compute neighbors and UMAP
    sc.pp.neighbors(adata, use_rep=rep_key)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)

    return adata, model

# Usage
adatas = {
    "study1": sc.read_h5ad("study1.h5ad"),
    "study2": sc.read_h5ad("study2.h5ad"),
    "study3": sc.read_h5ad("study3.h5ad")
}

adata_integrated, model = integrate_datasets(
    adatas,
    labels_key="cell_type"
)

sc.pl.umap(adata_integrated, color=["batch", "leiden", "cell_type"])
```

## 문제 해결

| 문제 | 원인 | 해결 방법 |
|------|------|-----------|
| 배치 간 혼합 부족 | 공유 유전자 부족 | HVG 수 증가, 유전자 중첩 확인 |
| 과도한 보정 | 생물학적 변이 제거 | scANVI에 레이블 적용 |
| 학습 발산 | 학습률 과다 | lr 감소, batch_size 증가 |
| NaN 손실 | 불량 데이터 | 전체 0인 세포/유전자 확인 |
| 메모리 오류 | 세포 수 과다 | batch_size 감소, GPU 사용 |
