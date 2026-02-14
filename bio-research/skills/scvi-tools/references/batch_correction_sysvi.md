# sysVI를 이용한 고급 배치 보정

이 레퍼런스는 주요 기술적 또는 연구 간 차이가 있는 데이터를 통합하기 위해 설계된 시스템 수준 배치 보정 방법인 sysVI를 다룹니다.

## 개요

sysVI(System Variational Inference)는 다음과 같은 시나리오에서 scVI를 확장합니다:
- 배치 효과가 매우 강한 경우 (서로 다른 기술)
- 표준 scVI가 생물학적 신호를 과도하게 보정하는 경우
- "시스템" 효과와 생물학적 변이를 분리해야 하는 경우

## sysVI와 scVI 사용 시점 비교

| 시나리오 | 권장 모델 |
|---------|----------|
| 동일 기술, 다른 샘플 | scVI |
| 10x v2 vs 10x v3 | scVI (일반적으로) |
| 10x vs Smart-seq2 | sysVI |
| 시퀀싱 깊이가 다른 경우 | 공변량이 포함된 scVI |
| 교차 연구 통합 | sysVI |
| 아틀라스 규모 통합 | sysVI |

## 사전 요구사항

```python
import scvi
import scanpy as sc
import numpy as np

print(f"scvi-tools version: {scvi.__version__}")
```

## sysVI 아키텍처 이해

sysVI는 변이를 다음과 같이 분리합니다:
1. **생물학적 변이**: 세포 유형, 상태, 궤적
2. **시스템 변이**: 기술, 연구, 실험실 효과

```
                    ┌─────────────────┐
Input counts ──────►│    Encoder      │
                    │                 │
System info ───────►│  (conditioned)  │
                    └────────┬────────┘
                             │
                    ┌────────▼────────┐
                    │   Latent z      │
                    │  (biological)   │
                    └────────┬────────┘
                             │
                    ┌────────▼────────┐
System info ───────►│    Decoder      │
                    │  (conditioned)  │
                    └────────┬────────┘
                             │
                    Reconstructed counts
```

## 기본 sysVI 워크플로우

### 1단계: 데이터 준비

```python
# Load datasets from different systems
adata1 = sc.read_h5ad("10x_data.h5ad")
adata2 = sc.read_h5ad("smartseq_data.h5ad")

# Add system labels
adata1.obs["system"] = "10x"
adata2.obs["system"] = "Smart-seq2"

# Add batch labels (within system)
# e.g., different samples within each technology

# Concatenate
adata = sc.concat([adata1, adata2])

# Store raw counts
adata.layers["counts"] = adata.X.copy()
```

### 2단계: HVG 선택

```python
# Select HVGs considering both batch and system
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=4000,  # More genes for cross-system
    flavor="seurat_v3",
    batch_key="system",  # Consider system for HVG
    layer="counts"
)

# Optionally: ensure overlap between systems
# Check HVGs are expressed in both systems
adata = adata[:, adata.var["highly_variable"]].copy()
```

### 3단계: sysVI 설정 및 학습

```python
# Setup AnnData
# Note: sysVI may be accessed differently depending on version
# Check scvi-tools documentation for current API

scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="sample",           # Within-system batches
    categorical_covariate_keys=["system"]  # System-level covariate
)

# For true sysVI (if available in your version)
# scvi.model.SysVI.setup_anndata(...)

# Create model with system awareness
model = scvi.model.SCVI(
    adata,
    n_latent=30,
    n_layers=2,
    gene_likelihood="nb"
)

# Train
model.train(max_epochs=300)
```

### 4단계: 표현 추출

```python
# Get latent representation
adata.obsm["X_integrated"] = model.get_latent_representation()

# Clustering and visualization
sc.pp.neighbors(adata, use_rep="X_integrated")
sc.tl.umap(adata)
sc.tl.leiden(adata)

# Check integration
sc.pl.umap(adata, color=["system", "leiden", "cell_type"])
```

## 대안: Harmony + scVI

교차 시스템 통합의 경우, 방법을 결합하면 좋은 결과를 얻을 수 있습니다:

```python
import scanpy.external as sce

# First run PCA
sc.pp.pca(adata)

# Apply Harmony for initial alignment
sce.pp.harmony_integrate(adata, key="system")

# Then train scVI on Harmony-corrected embedding
# Or use Harmony representation directly
```

## 대안: scVI에서 공변량 사용

중간 수준의 시스템 효과가 있는 경우:

```python
# Include system as categorical covariate
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="sample",
    categorical_covariate_keys=["system", "technology_version"]
)

model = scvi.model.SCVI(adata, n_latent=30)
model.train()
```

## 대안: 분리 모델 + 통합

매우 다른 시스템의 경우:

```python
# Train separate models
scvi.model.SCVI.setup_anndata(adata1, layer="counts", batch_key="sample")
model1 = scvi.model.SCVI(adata1)
model1.train()

scvi.model.SCVI.setup_anndata(adata2, layer="counts", batch_key="sample")
model2 = scvi.model.SCVI(adata2)
model2.train()

# Get latent spaces
adata1.obsm["X_scVI"] = model1.get_latent_representation()
adata2.obsm["X_scVI"] = model2.get_latent_representation()

# Align with CCA or Harmony
# ... additional alignment step
```

## 교차 시스템 통합 평가

### 시각적 평가

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Color by system
sc.pl.umap(adata, color="system", ax=axes[0], show=False, title="By System")

# Color by cell type
sc.pl.umap(adata, color="cell_type", ax=axes[1], show=False, title="By Cell Type")

# Color by expression of marker
sc.pl.umap(adata, color="CD3D", ax=axes[2], show=False, title="CD3D Expression")

plt.tight_layout()
```

### 정량적 메트릭

```python
# Using scib-metrics
from scib_metrics.benchmark import Benchmarker

bm = Benchmarker(
    adata,
    batch_key="system",
    label_key="cell_type",
    embedding_obsm_keys=["X_integrated"]
)

bm.benchmark()

# Key metrics:
# - Batch mixing (ASW_batch, Graph connectivity)
# - Bio conservation (NMI, ARI, ASW_label)
```

### LISI 점수

```python
# Local Inverse Simpson's Index
from scib_metrics import lisi

# Batch LISI (higher = better mixing)
batch_lisi = lisi.ilisi_graph(
    adata,
    batch_key="system",
    use_rep="X_integrated"
)

# Cell type LISI (lower = better preservation)
ct_lisi = lisi.clisi_graph(
    adata,
    label_key="cell_type",
    use_rep="X_integrated"
)

print(f"Batch LISI: {batch_lisi.mean():.3f}")
print(f"Cell type LISI: {ct_lisi.mean():.3f}")
```

## 특정 문제 상황 처리

### 유전자 세트가 다른 경우

```python
# Find common genes
common_genes = adata1.var_names.intersection(adata2.var_names)
print(f"Common genes: {len(common_genes)}")

# If too few, use gene mapping
# Or impute missing genes
```

### 시퀀싱 깊이가 다른 경우

```python
# Add depth as continuous covariate
adata.obs["log_counts"] = np.log1p(adata.obs["total_counts"])

scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="sample",
    continuous_covariate_keys=["log_counts"]
)
```

### 세포 유형 불균형

```python
# Check cell type distribution per system
import pandas as pd

ct_dist = pd.crosstab(adata.obs["system"], adata.obs["cell_type"], normalize="index")
print(ct_dist)

# If very unbalanced, consider:
# 1. Subsample to balance
# 2. Use scANVI with labels to preserve rare types
```

## 전체 파이프라인

```python
def integrate_cross_system(
    adatas: dict,
    system_key: str = "system",
    batch_key: str = "batch",
    cell_type_key: str = "cell_type",
    n_top_genes: int = 4000,
    n_latent: int = 30
):
    """
    Integrate datasets from different technological systems.

    Parameters
    ----------
    adatas : dict
        Dictionary of {system_name: AnnData}
    system_key : str
        Key for system annotation
    batch_key : str
        Key for within-system batch
    cell_type_key : str
        Key for cell type labels (optional)
    n_top_genes : int
        Number of HVGs
    n_latent : int
        Latent dimensions

    Returns
    -------
    Integrated AnnData with model
    """
    import scvi
    import scanpy as sc

    # Add system labels and concatenate
    for system_name, adata in adatas.items():
        adata.obs[system_key] = system_name

    adata = sc.concat(list(adatas.values()))

    # Find common genes
    for name, ad in adatas.items():
        if name == list(adatas.keys())[0]:
            common_genes = set(ad.var_names)
        else:
            common_genes = common_genes.intersection(ad.var_names)

    adata = adata[:, list(common_genes)].copy()
    print(f"Common genes: {len(common_genes)}")

    # Store counts
    adata.layers["counts"] = adata.X.copy()

    # HVG selection
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor="seurat_v3",
        batch_key=system_key,
        layer="counts"
    )
    adata = adata[:, adata.var["highly_variable"]].copy()

    # Setup with system as covariate
    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key=batch_key if batch_key in adata.obs else None,
        categorical_covariate_keys=[system_key]
    )

    # Train
    model = scvi.model.SCVI(adata, n_latent=n_latent, n_layers=2)
    model.train(max_epochs=300, early_stopping=True)

    # Get representation
    adata.obsm["X_integrated"] = model.get_latent_representation()

    # Clustering
    sc.pp.neighbors(adata, use_rep="X_integrated")
    sc.tl.umap(adata)
    sc.tl.leiden(adata)

    return adata, model

# Usage
adatas = {
    "10x_v3": sc.read_h5ad("10x_v3_data.h5ad"),
    "Smart-seq2": sc.read_h5ad("smartseq_data.h5ad"),
    "Drop-seq": sc.read_h5ad("dropseq_data.h5ad")
}

adata_integrated, model = integrate_cross_system(adatas)

# Visualize
sc.pl.umap(adata_integrated, color=["system", "leiden"])
```

## 문제 해결

| 문제 | 원인 | 해결 방법 |
|------|------|----------|
| 시스템 간 혼합 안 됨 | 효과가 너무 강함 | 유전자 수 늘리기, n_latent 증가 |
| 과도한 보정 | 모델이 너무 공격적 | n_layers 줄이기, scANVI 사용 |
| 공통 유전자 부족 | 다른 플랫폼 | 유전자 이름 매핑 사용 |
| 하나의 시스템이 지배적 | 크기 불균형 | 큰 데이터셋 서브샘플링 |

## 주요 참고문헌

- Lopez et al. (2018) "Deep generative modeling for single-cell transcriptomics"
- Luecken et al. (2022) "Benchmarking atlas-level data integration in single-cell genomics"
