# scvi-tools를 위한 데이터 준비

이 레퍼런스는 scvi-tools 모델에서 사용할 AnnData 객체를 올바르게 준비하는 방법을 다룹니다.

## 개요

scvi-tools에서 적절한 데이터 준비는 매우 중요합니다. 핵심 요구사항은 다음과 같습니다:
1. **원시 카운트** (정규화되지 않은 상태)
2. **고변이 유전자 선택**
3. **적절한 setup_anndata() 호출**

## 1단계: 데이터 로드 및 검사

```python
import scanpy as sc
import scvi
import numpy as np

# Load data
adata = sc.read_h5ad("data.h5ad")

# Check what's in adata.X
print(f"Shape: {adata.shape}")
print(f"X dtype: {adata.X.dtype}")
print(f"X contains integers: {np.allclose(adata.X.data, adata.X.data.astype(int))}")
print(f"X min: {adata.X.min()}, max: {adata.X.max()}")
```

### 원시 카운트 확인

```python
# scvi-tools needs INTEGER counts
# If X appears normalized, check for raw counts

if hasattr(adata, 'raw') and adata.raw is not None:
    print("Found adata.raw")
    # Use raw counts
    adata = adata.raw.to_adata()

# Or check layers
if 'counts' in adata.layers:
    print("Found counts layer")
    # Will specify layer in setup_anndata
```

## 2단계: 기본 필터링

```python
# Filter cells (standard QC)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, max_genes=5000)

# Calculate mito percent if not present
# Handle both human (MT-) and mouse (mt-, Mt-) mitochondrial genes
adata.var['mt'] = (
    adata.var_names.str.startswith('MT-') |
    adata.var_names.str.startswith('mt-') |
    adata.var_names.str.startswith('Mt-')
)
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
adata = adata[adata.obs['pct_counts_mt'] < 20].copy()

# Filter genes
sc.pp.filter_genes(adata, min_cells=3)

print(f"After filtering: {adata.shape}")
```

## 3단계: 원시 카운트 저장

**중요**: 정규화 전에 반드시 원시 카운트를 보존해야 합니다.

```python
# Store raw counts in a layer
adata.layers["counts"] = adata.X.copy()

# Now you can normalize for other purposes (HVG selection)
# But scvi will use the counts layer
```

## 4단계: 고변이 유전자(HVG) 선택

scvi-tools는 1,500~5,000개의 HVG에서 최적의 성능을 보입니다.

### 단일 배치 데이터

```python
# Normalize for HVG selection only
adata_hvg = adata.copy()
sc.pp.normalize_total(adata_hvg, target_sum=1e4)
sc.pp.log1p(adata_hvg)

# Select HVGs
sc.pp.highly_variable_genes(
    adata_hvg,
    n_top_genes=2000,
    flavor="seurat"  # or "cell_ranger"
)

# Transfer HVG annotation
adata.var['highly_variable'] = adata_hvg.var['highly_variable']
```

### 다중 배치 데이터 (권장)

```python
# Use seurat_v3 flavor with batch_key
# This selects genes variable across batches
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    flavor="seurat_v3",
    batch_key="batch",  # Your batch column
    layer="counts"      # Use raw counts
)
```

### HVG로 서브셋

```python
# Subset to highly variable genes
adata = adata[:, adata.var['highly_variable']].copy()
print(f"After HVG selection: {adata.shape}")
```

## 5단계: AnnData 설정

`setup_anndata()` 함수는 모델에 사용할 데이터를 등록합니다.

### 기본 설정

```python
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts"  # Specify layer with raw counts
)
```

### 배치 정보 포함

```python
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="batch"  # Column in adata.obs
)
```

### 세포 유형 라벨 포함 (scANVI용)

```python
scvi.model.SCANVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="batch",
    labels_key="cell_type"  # Column with cell type labels
)
```

### 연속형 공변량 포함

```python
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="batch",
    continuous_covariate_keys=["percent_mito", "n_genes"]
)
```

### 범주형 공변량 포함

```python
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="batch",
    categorical_covariate_keys=["donor", "technology"]
)
```

## 다중 모달 데이터 설정

### CITE-seq (totalVI용)

```python
# Protein data in adata.obsm
# RNA in adata.X, protein in separate matrix

# Add protein data
adata.obsm["protein_expression"] = protein_counts  # numpy array

# Setup for totalVI
scvi.model.TOTALVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="batch",
    protein_expression_obsm_key="protein_expression"
)
```

### 멀티옴 RNA+ATAC (MultiVI용)

```python
# RNA and ATAC in separate AnnData objects or MuData

import mudata as md

# If using MuData
mdata = md.read("multiome.h5mu")

scvi.model.MULTIVI.setup_mudata(
    mdata,
    rna_layer="counts",
    protein_layer=None,
    batch_key="batch",
    modalities={"rna": "rna", "accessibility": "atac"}
)
```

## 전체 준비 파이프라인

전체 준비 함수는 `scripts/model_utils.py`의 `prepare_adata()`를 사용합니다:

```python
from model_utils import prepare_adata

# Prepare data with QC, HVG selection, and layer setup
adata = prepare_adata(
    adata,
    batch_key="batch",
    n_top_genes=2000,
    min_genes=200,
    max_mito_pct=20
)

# Then setup for your model
import scvi
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
```

이 함수는 다음을 처리합니다:
- 미토콘드리아 QC 필터링
- 세포 및 유전자 필터링
- 레이어에 카운트 저장
- HVG 선택 (batch_key가 제공된 경우 배치 인식)
- HVG로 서브셋

## 설정 확인

```python
# View registered data
print(adata.uns['_scvi_manager_uuid'])
print(adata.uns['_scvi_adata_minify_type'])

# For scVI
scvi.model.SCVI.view_anndata_setup(adata)
```

## 일반적인 문제 및 해결 방법

| 문제 | 원인 | 해결 방법 |
|------|------|----------|
| "X should contain integers" | X에 정규화된 데이터가 있음 | layer="counts" 사용 |
| "batch_key not found" | 잘못된 열 이름 | adata.obs.columns 확인 |
| 희소 행렬 오류 | 호환되지 않는 형식 | 변환: adata.X = adata.X.toarray() |
| 메모리 오류 | 유전자 수 과다 | 먼저 HVG로 서브셋 |
| 데이터에 NaN 존재 | 결측값 | 필터링 또는 결측값 대체 |

## 데이터 형식 참조

### 필수 항목

- `adata.X` 또는 `adata.layers["counts"]`: 원시 정수 카운트 (희소 형식 가능)
- `adata.obs`: 세포 메타데이터 DataFrame
- `adata.var`: 유전자 메타데이터 DataFrame

### 권장 항목

- `adata.obs["batch"]`: 배치/샘플 식별자
- `adata.var["highly_variable"]`: HVG 불리언 마스크

### scANVI용

- `adata.obs["labels"]`: 세포 유형 주석
- 라벨이 없는 세포에 대해 "Unknown"을 포함할 수 있습니다
