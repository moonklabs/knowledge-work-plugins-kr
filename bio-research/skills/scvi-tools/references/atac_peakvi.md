# PeakVI를 이용한 scATAC-seq 분석

이 레퍼런스는 PeakVI를 사용한 단일세포 ATAC-seq 분석의 차원 축소, 배치 보정 및 차등 접근성 분석을 다룹니다.

## 개요

PeakVI는 scATAC-seq 데이터를 위한 딥 생성 모델로 다음과 같은 기능을 제공합니다:
- 이진 접근성(피크 열림/닫힘) 모델링
- 배치 효과 처리
- 클러스터링을 위한 잠재 표현 제공
- 차등 접근성 분석 수행

## 사전 요구사항

```python
import scvi
import scanpy as sc
import numpy as np
import anndata as ad

print(f"scvi-tools version: {scvi.__version__}")
```

## 1단계: ATAC 데이터 로드 및 준비

### 10x Genomics (Cell Ranger ATAC)에서 로드

```python
# Peak-cell matrix from fragments
# Usually in filtered_peak_bc_matrix format

adata = sc.read_10x_h5("filtered_peak_bc_matrix.h5")

# Or from mtx format
adata = sc.read_10x_mtx("filtered_peak_bc_matrix/")

# Check structure
print(f"Cells: {adata.n_obs}, Peaks: {adata.n_vars}")
print(f"Sparsity: {1 - adata.X.nnz / (adata.n_obs * adata.n_vars):.2%}")
```

### ArchR/Signac에서 로드

```python
# Export from ArchR (in R)
# saveArchRProject(proj, outputDirectory="atac_export", load=FALSE)
# Then read the exported files in Python

# From Signac:
# Export peak matrix and metadata
```

## 2단계: 품질 관리

```python
# Calculate QC metrics
sc.pp.calculate_qc_metrics(adata, inplace=True)

# Key metrics for ATAC:
# - n_genes_by_counts: peaks per cell (should rename)
# - total_counts: fragments per cell
adata.obs['n_peaks'] = adata.obs['n_genes_by_counts']
adata.obs['total_fragments'] = adata.obs['total_counts']

# Filter cells
adata = adata[adata.obs['n_peaks'] > 500].copy()
adata = adata[adata.obs['n_peaks'] < 50000].copy()  # Remove potential doublets

# Filter peaks (accessible in at least n cells)
sc.pp.filter_genes(adata, min_cells=10)

print(f"After QC: {adata.shape}")
```

### 데이터 이진화

```python
# PeakVI works with binary accessibility
# Binarize if not already binary
adata.X = (adata.X > 0).astype(np.float32)

# Verify
print(f"Unique values: {np.unique(adata.X.data)}")
```

## 3단계: 피처 선택

RNA-seq와 달리 ATAC에서의 피크 선택은 아직 표준화되지 않았습니다. 다음과 같은 옵션이 있습니다:

### 옵션 A: 가장 접근 가능한 피크

```python
# Select top peaks by accessibility frequency
peak_accessibility = np.array(adata.X.sum(axis=0)).flatten()
top_peaks = np.argsort(peak_accessibility)[-50000:]  # Top 50k peaks

adata = adata[:, top_peaks].copy()
```

### 옵션 B: 가변 피크

```python
# Select peaks with high variance
# (Most informative for clustering)
from sklearn.feature_selection import VarianceThreshold

selector = VarianceThreshold(threshold=0.05)
selector.fit(adata.X)
adata = adata[:, selector.get_support()].copy()
```

### 옵션 C: 유전자 근처 피크

```python
# Keep peaks within promoter regions or gene bodies
# Requires peak annotation
# gene_peaks = peaks with gene annotation
# adata = adata[:, adata.var['near_gene']].copy()
```

## 4단계: 배치 정보 추가

```python
# Add batch annotation if multiple samples
adata.obs['batch'] = adata.obs['sample_id']  # Or appropriate column

print(adata.obs['batch'].value_counts())
```

## 5단계: PeakVI 설정 및 학습

```python
# Setup AnnData
scvi.model.PEAKVI.setup_anndata(
    adata,
    batch_key="batch"  # Optional, omit for single batch
)

# Create model
model = scvi.model.PEAKVI(
    adata,
    n_latent=20,      # Latent dimensions
    n_layers_encoder=2,
    n_layers_decoder=2
)

# Train
model.train(
    max_epochs=200,
    early_stopping=True,
    batch_size=128
)

# Check training
model.history['elbo_train'].plot()
```

## 6단계: 잠재 표현 획득

```python
# Latent space for downstream analysis
adata.obsm["X_PeakVI"] = model.get_latent_representation()

# Clustering and visualization
sc.pp.neighbors(adata, use_rep="X_PeakVI", n_neighbors=15)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)

# Visualize
sc.pl.umap(adata, color=['leiden', 'batch'], ncols=2)
```

## 7단계: 차등 접근성 분석

```python
# Differential accessibility between clusters
da_results = model.differential_accessibility(
    groupby='leiden',
    group1='0',
    group2='1'
)

# Filter significant peaks
da_sig = da_results[
    (da_results['is_da_fdr_0.05']) &
    (abs(da_results['lfc_mean']) > 1)
]

print(f"Significant DA peaks: {len(da_sig)}")
print(da_sig.head())
```

### 조건 간 차등 접근성

```python
# Compare conditions within cell type
adata_subset = adata[adata.obs['cell_type'] == 'CD4 T cells'].copy()

da_condition = model.differential_accessibility(
    groupby='condition',
    group1='treated',
    group2='control'
)
```

## 8단계: 피크 주석

```python
# Annotate peaks with nearest genes
# Using pybedtools or similar

# Example peak name format: chr1:1000-2000
# Parse into bed format for annotation

import pandas as pd

def parse_peak_names(peak_names):
    """Parse peak names into bed format."""
    records = []
    for peak in peak_names:
        chrom, coords = peak.split(':')
        start, end = coords.split('-')
        records.append({
            'chrom': chrom,
            'start': int(start),
            'end': int(end),
            'peak': peak
        })
    return pd.DataFrame(records)

peak_bed = parse_peak_names(adata.var_names)
```

## 9단계: 모티프 분석

```python
# Export significant peaks for motif analysis
# Use HOMER, MEME, or chromVAR

# Export peak sequences
sig_peaks = da_sig.index.tolist()
peak_bed_sig = peak_bed[peak_bed['peak'].isin(sig_peaks)]
peak_bed_sig.to_csv("significant_peaks.bed", sep='\t', index=False, header=False)

# Then run HOMER:
# findMotifsGenome.pl significant_peaks.bed hg38 motif_output/ -size 200
```

## 10단계: 유전자 활성도 점수

```python
# Compute gene activity from peak accessibility
# (Requires peak-gene annotations)

def compute_gene_activity(adata, peak_gene_map):
    """
    Compute gene activity scores from peak accessibility.

    Parameters
    ----------
    adata : AnnData
        ATAC data with peaks
    peak_gene_map : dict
        Mapping of peaks to genes

    Returns
    -------
    AnnData with gene activity scores
    """
    from scipy.sparse import csr_matrix

    genes = list(set(peak_gene_map.values()))
    gene_matrix = np.zeros((adata.n_obs, len(genes)))

    for i, gene in enumerate(genes):
        gene_peaks = [p for p, g in peak_gene_map.items() if g == gene]
        if gene_peaks:
            peak_idx = [list(adata.var_names).index(p) for p in gene_peaks if p in adata.var_names]
            if peak_idx:
                gene_matrix[:, i] = np.array(adata.X[:, peak_idx].sum(axis=1)).flatten()

    adata_gene = ad.AnnData(
        X=csr_matrix(gene_matrix),
        obs=adata.obs.copy(),
        var=pd.DataFrame(index=genes)
    )

    return adata_gene
```

## 전체 파이프라인

```python
def analyze_scatac(
    adata,
    batch_key=None,
    n_top_peaks=50000,
    n_latent=20,
    resolution=0.5
):
    """
    Complete scATAC-seq analysis with PeakVI.

    Parameters
    ----------
    adata : AnnData
        Raw peak-cell matrix
    batch_key : str, optional
        Batch annotation column
    n_top_peaks : int
        Number of top peaks to use
    n_latent : int
        Latent dimensions
    resolution : float
        Leiden clustering resolution

    Returns
    -------
    Tuple of (processed AnnData, trained model)
    """
    import scvi
    import scanpy as sc
    import numpy as np

    adata = adata.copy()

    # QC
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata = adata[adata.obs['n_genes_by_counts'] > 500].copy()
    sc.pp.filter_genes(adata, min_cells=10)

    # Binarize
    adata.X = (adata.X > 0).astype(np.float32)

    # Select top peaks
    if adata.n_vars > n_top_peaks:
        peak_accessibility = np.array(adata.X.sum(axis=0)).flatten()
        top_peaks = np.argsort(peak_accessibility)[-n_top_peaks:]
        adata = adata[:, top_peaks].copy()

    # Setup PeakVI
    scvi.model.PEAKVI.setup_anndata(adata, batch_key=batch_key)

    # Train
    model = scvi.model.PEAKVI(adata, n_latent=n_latent)
    model.train(max_epochs=200, early_stopping=True)

    # Latent representation
    adata.obsm["X_PeakVI"] = model.get_latent_representation()

    # Clustering
    sc.pp.neighbors(adata, use_rep="X_PeakVI")
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=resolution)

    return adata, model

# Usage
adata, model = analyze_scatac(
    adata,
    batch_key="sample",
    n_top_peaks=50000
)

# Visualize
sc.pl.umap(adata, color=['leiden', 'sample'])

# Differential accessibility
da_results = model.differential_accessibility(
    groupby='leiden',
    group1='0',
    group2='1'
)
```

## scRNA-seq와의 통합

멀티옴 데이터 또는 동일 세포에서 분리된 RNA/ATAC 데이터의 경우:

```python
# See MultiVI for joint RNA+ATAC analysis
# Or use WNN (weighted nearest neighbors) approach

# Transfer labels from RNA to ATAC using shared latent space
```

## 문제 해결

| 문제 | 원인 | 해결 방법 |
|------|------|----------|
| 학습 속도 느림 | 피크 수가 너무 많음 | 상위 50k 피크로 서브셋 |
| 클러스터링 품질 저하 | 정보량이 높은 피크 부족 | 가변 피크 사용 |
| 배치가 지배적 | 강한 기술적 효과 | batch_key 설정 확인 |
| 메모리 오류 | 피크 행렬 크기가 큼 | 희소 형식 사용, 피크 수 줄이기 |

## 주요 참고문헌

- Ashuach et al. (2022) "PeakVI: A deep generative model for single-cell chromatin accessibility analysis"
