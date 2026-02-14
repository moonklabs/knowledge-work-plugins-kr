# scvi-tools 문제 해결 가이드

본 레퍼런스는 모든 scvi-tools 모델에 걸친 일반적인 문제의 진단 및 해결을 위한 통합 가이드를 제공합니다.

## 빠른 진단

| 증상 | 예상 원인 | 빠른 해결 |
|------|----------|-----------|
| "X should contain integers" | X에 정규화된 데이터가 있음 | setup에서 `layer="counts"` 사용 |
| CUDA out of memory | GPU 메모리 부족 | `batch_size` 감소, 더 작은 모델 사용 |
| 학습 손실이 NaN | 불량 데이터 또는 학습률 문제 | 전체 0인 세포/유전자 확인 |
| 배치 간 혼합 부족 | 공유 피처 부족 | HVG 수 증가, 유전자 중첩 확인 |
| 과도한 보정 | 너무 공격적인 통합 | scANVI에 레이블 적용 |
| Import 오류 | 의존성 누락 | `pip install scvi-tools[all]` |

## 데이터 형식 문제

### 문제: Seurat에서 가져온 CITE-seq 단백질 데이터가 CLR 정규화됨

**원인**: Seurat의 `NormalizeData(normalization.method = "CLR")`이 원시 ADT 카운트를 변환합니다. totalVI는 단백질 데이터에 원시 정수 카운트가 필요합니다.

**증상**:
- 단백질 값이 정수가 아님
- 단백질 값에 음수 포함
- 모델 학습 결과가 저조

**해결 방법**:
```python
# Check if protein data is normalized
protein = adata.obsm["protein_expression"]
print(f"Min value: {protein.min()}")  # Should be 0 if raw counts
print(f"Contains integers: {np.allclose(protein, protein.astype(int))}")

# If importing from Seurat, use the raw counts assay, not the normalized one
# In R/Seurat, export the RNA assay's counts slot, not the data slot
# GetAssayData(seurat_obj, assay = "ADT", slot = "counts")
```

### 문제: "layer not found" 또는 "X should contain integers"

**원인**: scvi-tools는 정규화된 데이터가 아닌 원시 정수 카운트가 필요합니다.

**해결 방법**:
```python
# Check if X contains integers
import numpy as np
print(f"X max: {adata.X.max()}")
print(f"Contains integers: {np.allclose(adata.X.data, adata.X.data.astype(int))}")

# If normalized, recover from raw
if hasattr(adata, 'raw') and adata.raw is not None:
    adata = adata.raw.to_adata()

# Or use existing counts layer
adata.layers["counts"] = adata.X.copy()
scvi.model.SCVI.setup_anndata(adata, layer="counts")
```

### 문제: 희소 행렬 오류

**원인**: 호환되지 않는 희소 형식이거나 밀집 배열이 필요합니다.

**해결 방법**:
```python
from scipy.sparse import csr_matrix

# Convert to CSR format (most compatible)
if hasattr(adata.X, 'toarray'):
    adata.X = csr_matrix(adata.X)

# Or convert to dense if small enough
if adata.n_obs * adata.n_vars < 1e8:
    adata.X = adata.X.toarray()
```

### 문제: 데이터에 NaN 또는 Inf 값 존재

**원인**: 결측값 또는 손상된 데이터입니다.

**해결 방법**:
```python
import numpy as np

# Check for issues
X = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
print(f"NaN count: {np.isnan(X).sum()}")
print(f"Inf count: {np.isinf(X).sum()}")
print(f"Negative count: {(X < 0).sum()}")

# Replace NaN/Inf with 0
X = np.nan_to_num(X, nan=0, posinf=0, neginf=0)
X = np.clip(X, 0, None)  # Ensure non-negative
adata.X = csr_matrix(X)
```

### 문제: batch_key 또는 labels_key를 찾을 수 없음

**원인**: adata.obs의 열 이름 불일치입니다.

**해결 방법**:
```python
# List available columns
print(adata.obs.columns.tolist())

# Check for similar names
for col in adata.obs.columns:
    if 'batch' in col.lower() or 'sample' in col.lower():
        print(f"Potential batch column: {col}")
```

## GPU 및 메모리 문제

### 문제: CUDA out of memory

**원인**: 모델 또는 배치가 GPU 메모리에 맞지 않습니다.

**해결 방법** (순서대로 시도):

```python
# 1. Reduce batch size
model.train(batch_size=64)  # Default is 128

# 2. Use smaller model architecture
model = scvi.model.SCVI(
    adata,
    n_latent=10,   # Default is 10-30
    n_layers=1     # Default is 1-2
)

# 3. Subset to fewer genes
sc.pp.highly_variable_genes(adata, n_top_genes=1500)
adata = adata[:, adata.var['highly_variable']].copy()

# 4. Clear GPU cache between models
import torch
torch.cuda.empty_cache()

# 5. Use CPU if GPU is too small
model.train(accelerator="cpu")
```

### 문제: GPU가 감지되지 않음

**원인**: CUDA가 설치되지 않았거나 버전이 불일치합니다.

**진단**:
```python
import torch
print(f"CUDA available: {torch.cuda.is_available()}")
print(f"PyTorch version: {torch.__version__}")
print(f"CUDA version: {torch.version.cuda}")
```

**해결 방법**:
```bash
# Check system CUDA
nvidia-smi
nvcc --version

# Reinstall PyTorch with matching CUDA
pip install torch --index-url https://download.pytorch.org/whl/cu118  # For CUDA 11.8
# Or
pip install torch --index-url https://download.pytorch.org/whl/cu121  # For CUDA 12.1
```

### 문제: 대규모 데이터셋에서 메모리 오류

**원인**: 데이터셋이 시스템 RAM에 비해 너무 큽니다.

**해결 방법**:
```python
# 1. Process in chunks (for very large data)
# Subsample for initial exploration
adata_sample = adata[np.random.choice(adata.n_obs, 50000, replace=False)].copy()

# 2. Use backed mode for AnnData
adata = sc.read_h5ad("large_data.h5ad", backed='r')

# 3. Reduce gene count aggressively
adata = adata[:, adata.var['highly_variable']].copy()
```

## 학습 문제

### 문제: 학습 손실이 NaN

**원인**: 수치 불안정, 불량 데이터 또는 학습률 문제입니다.

**해결 방법**:
```python
# 1. Check for problematic cells/genes
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# 2. Remove cells with zero counts
adata = adata[adata.X.sum(axis=1) > 0].copy()

# 3. Use gradient clipping (built into scvi-tools)
model.train(max_epochs=200, early_stopping=True)
```

### 문제: 학습이 수렴하지 않음

**원인**: 에폭 부족, 부적절한 하이퍼파라미터 또는 데이터 문제입니다.

**해결 방법**:
```python
# 1. Train longer
model.train(max_epochs=400)

# 2. Check training curves
import matplotlib.pyplot as plt
plt.plot(model.history['elbo_train'])
plt.plot(model.history['elbo_validation'])
plt.xlabel('Epoch')
plt.ylabel('ELBO')
plt.legend(['Train', 'Validation'])

# 3. Adjust model size for data size
# Small data (<10k cells): smaller model
model = scvi.model.SCVI(adata, n_latent=10, n_layers=1, dropout_rate=0.2)

# Large data (>100k cells): can use larger model
model = scvi.model.SCVI(adata, n_latent=30, n_layers=2)
```

### 문제: 과적합 (검증 손실 증가)

**원인**: 모델이 너무 복잡하거나 학습 기간이 너무 깁니다.

**해결 방법**:
```python
# 1. Enable early stopping
model.train(early_stopping=True, early_stopping_patience=10)

# 2. Add regularization
model = scvi.model.SCVI(adata, dropout_rate=0.2)

# 3. Reduce model complexity
model = scvi.model.SCVI(adata, n_layers=1)
```

## 통합 문제

### 문제: 배치 간 혼합 부족

**원인**: 공유 피처 부족, 강한 생물학적 차이 또는 기술적 문제입니다.

**해결 방법**:
```python
# 1. Check gene overlap between batches
for batch in adata.obs['batch'].unique():
    batch_genes = adata[adata.obs['batch'] == batch].var_names
    print(f"{batch}: {len(batch_genes)} genes")

# 2. Use more HVGs
sc.pp.highly_variable_genes(adata, n_top_genes=4000, batch_key="batch")

# 3. Train longer
model.train(max_epochs=400)

# 4. Increase latent dimensions
model = scvi.model.SCVI(adata, n_latent=50)
```

### 문제: 과도한 보정 (생물학적 신호 손실)

**원인**: 모델이 너무 많은 변이를 제거합니다.

**해결 방법**:
```python
# 1. Use scANVI with cell type labels
scvi.model.SCANVI.from_scvi_model(scvi_model, labels_key="cell_type")

# 2. Reduce model capacity
model = scvi.model.SCVI(adata, n_latent=10)

# 3. Use categorical covariates instead of batch_key
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    categorical_covariate_keys=["batch"]  # Less aggressive than batch_key
)
```

### 문제: 하나의 배치가 클러스터를 지배

**원인**: 배치 크기 불균형 또는 불완전한 통합입니다.

**해결 방법**:
```python
# 1. Check batch distribution
print(adata.obs['batch'].value_counts())

# 2. Subsample to balance
from sklearn.utils import resample
balanced = []
min_size = adata.obs['batch'].value_counts().min()
for batch in adata.obs['batch'].unique():
    batch_data = adata[adata.obs['batch'] == batch]
    balanced.append(batch_data[np.random.choice(len(batch_data), min_size, replace=False)])
adata_balanced = sc.concat(balanced)
```

## 모델별 문제

### scANVI: 레이블 전이 품질 저하

**해결 방법**:
```python
# 1. Check label distribution
print(adata.obs['cell_type'].value_counts())

# 2. Use Unknown for low-confidence cells
adata.obs.loc[adata.obs['prediction_score'] < 0.5, 'cell_type'] = 'Unknown'

# 3. Train scVI longer before scANVI
scvi_model.train(max_epochs=300)
scanvi_model = scvi.model.SCANVI.from_scvi_model(scvi_model, labels_key="cell_type")
scanvi_model.train(max_epochs=100)
```

### totalVI: 노이즈가 많은 단백질 신호

**해결 방법**:
```python
# 1. Use denoised protein values
_, protein_denoised = model.get_normalized_expression(return_mean=True)

# 2. Check isotype controls
# Isotype controls should have low expression
for i, name in enumerate(adata.uns["protein_names"]):
    if 'isotype' in name.lower():
        print(f"{name}: mean={adata.obsm['protein_expression'][:, i].mean():.1f}")
```

### PeakVI: 클러스터링 품질 저하

**해결 방법**:
```python
# 1. Use more variable peaks
from sklearn.feature_selection import VarianceThreshold
selector = VarianceThreshold(threshold=0.05)
adata = adata[:, selector.fit(adata.X).get_support()].copy()

# 2. Binarize data
adata.X = (adata.X > 0).astype(np.float32)
```

### MultiVI: 모달리티 간 세포 수 불일치

**해결 방법**:
```python
# Ensure same cells in same order
common_cells = adata_rna.obs_names.intersection(adata_atac.obs_names)
adata_rna = adata_rna[common_cells].copy()
adata_atac = adata_atac[common_cells].copy()
```

### DestVI: 디컨볼루션 품질 저하

**해결 방법**:
```python
# 1. Check gene overlap
common_genes = adata_ref.var_names.intersection(adata_spatial.var_names)
print(f"Common genes: {len(common_genes)}")  # Should be >1000

# 2. Use tissue-matched reference
# Reference should contain all cell types expected in spatial data

# 3. Check reference quality
print(adata_ref.obs['cell_type'].value_counts())
```

## 버전 호환성

### scvi-tools 1.x vs 0.x API 변경 사항

주요 차이점:
```python
# 0.x API
scvi.data.setup_anndata(adata, ...)

# 1.x API (current)
scvi.model.SCVI.setup_anndata(adata, ...)
```

### 버전 확인
```python
import scvi
import scanpy as sc
import anndata
import torch

print(f"scvi-tools: {scvi.__version__}")
print(f"scanpy: {sc.__version__}")
print(f"anndata: {anndata.__version__}")
print(f"torch: {torch.__version__}")
```

### 권장 버전 (2024년 후반 기준)
```
scvi-tools>=1.0.0
scanpy>=1.9.0
anndata>=0.9.0
torch>=2.0.0
```

## 도움 받기

1. **문서 확인**: https://docs.scvi-tools.org/
2. **GitHub 이슈**: https://github.com/scverse/scvi-tools/issues
3. **Discourse 포럼**: https://discourse.scverse.org/
4. **튜토리얼**: https://docs.scvi-tools.org/en/stable/tutorials/index.html

이슈 보고 시 포함할 내용:
- scvi-tools 버전 (`scvi.__version__`)
- Python 버전
- 전체 오류 트레이스백
- 최소한의 재현 가능 예제
