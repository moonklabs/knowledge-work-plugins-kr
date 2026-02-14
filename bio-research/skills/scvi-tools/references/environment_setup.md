# scvi-tools 환경 설정

이 레퍼런스는 scvi-tools의 설치 및 환경 구성을 다룹니다.

## 설치 옵션

### 옵션 1: Conda 환경 (권장)

```bash
# Create environment with GPU support
conda create -n scvi-env python=3.10
conda activate scvi-env

# Install scvi-tools
pip install scvi-tools

# For GPU acceleration (recommended for large datasets)
pip install torch --index-url https://download.pytorch.org/whl/cu118

# Common dependencies
pip install scanpy leidenalg
```

### 옵션 2: Pip만 사용

```bash
# Create virtual environment
python -m venv scvi-env
source scvi-env/bin/activate  # Linux/Mac
# scvi-env\Scripts\activate   # Windows

# Install
pip install scvi-tools scanpy
```

### 옵션 3: 공간 분석 지원 포함

```bash
conda create -n scvi-spatial python=3.10
conda activate scvi-spatial

pip install scvi-tools scanpy squidpy
```

### 옵션 4: MuData 지원 포함 (멀티옴)

```bash
pip install scvi-tools mudata muon
```

## 설치 확인

```python
import scvi
import torch
import scanpy as sc

print(f"scvi-tools version: {scvi.__version__}")
print(f"scanpy version: {sc.__version__}")
print(f"PyTorch version: {torch.__version__}")
print(f"GPU available: {torch.cuda.is_available()}")

if torch.cuda.is_available():
    print(f"GPU device: {torch.cuda.get_device_name(0)}")
    print(f"GPU memory: {torch.cuda.get_device_properties(0).total_mem / 1e9:.1f} GB")
```

## GPU 구성

### CUDA 버전 확인

```bash
nvidia-smi
nvcc --version
```

### PyTorch CUDA 버전

| CUDA 버전 | PyTorch 설치 명령어 |
|-----------|-------------------|
| CUDA 11.8 | `pip install torch --index-url https://download.pytorch.org/whl/cu118` |
| CUDA 12.1 | `pip install torch --index-url https://download.pytorch.org/whl/cu121` |
| CPU만 사용 | `pip install torch --index-url https://download.pytorch.org/whl/cpu` |

### 메모리 관리

```python
import torch

# Clear GPU cache between models
torch.cuda.empty_cache()

# Monitor memory usage
print(f"Allocated: {torch.cuda.memory_allocated() / 1e9:.2f} GB")
print(f"Cached: {torch.cuda.memory_reserved() / 1e9:.2f} GB")
```

## 일반적인 문제

| 문제 | 원인 | 해결 방법 |
|------|------|----------|
| `CUDA out of memory` | GPU 메모리 부족 | batch_size 줄이기, 더 작은 모델 사용 |
| `No GPU detected` | CUDA 미설치 | PyTorch에 맞는 CUDA 툴킷 설치 |
| `Version mismatch` | PyTorch/CUDA 비호환 | 올바른 CUDA 버전의 PyTorch 재설치 |
| `Import error scvi` | 의존성 누락 | `pip install scvi-tools[all]` |

## Jupyter 설정

```bash
# Install Jupyter kernel
pip install ipykernel
python -m ipykernel install --user --name scvi-env --display-name "scvi-tools"

# For interactive plots
pip install matplotlib seaborn
```

## 권장 패키지 버전

재현성을 위해 버전을 고정합니다:

```bash
pip install \
    scvi-tools>=1.0.0 \
    scanpy>=1.9.0 \
    anndata>=0.9.0 \
    torch>=2.0.0
```

## 버전 호환성 가이드

### scvi-tools 1.x vs 0.x API 변경사항

1.x 릴리스에서 호환성을 깨는 변경사항이 도입되었습니다. 주요 차이점은 다음과 같습니다:

| 작업 | 0.x API (더 이상 사용되지 않음) | 1.x API (현재) |
|-----|-------------------------------|----------------|
| 데이터 설정 | `scvi.data.setup_anndata(adata, ...)` | `scvi.model.SCVI.setup_anndata(adata, ...)` |
| 데이터 등록 | `scvi.data.register_tensor_from_anndata(...)` | `setup_anndata`에 내장 |
| 설정 보기 | `scvi.data.view_anndata_setup(adata)` | `scvi.model.SCVI.view_anndata_setup(adata)` |

### 0.x에서 1.x로 마이그레이션

```python
# OLD (0.x) - DEPRECATED
import scvi
scvi.data.setup_anndata(adata, layer="counts", batch_key="batch")
model = scvi.model.SCVI(adata)

# NEW (1.x) - CURRENT
import scvi
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
model = scvi.model.SCVI(adata)
```

### 모델별 설정 (1.x)

각 모델에는 고유한 설정 메서드가 있습니다:

```python
# scVI
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")

# scANVI
scvi.model.SCANVI.setup_anndata(adata, layer="counts", batch_key="batch", labels_key="cell_type")

# totalVI
scvi.model.TOTALVI.setup_anndata(adata, layer="counts", protein_expression_obsm_key="protein")

# MultiVI (uses MuData)
scvi.model.MULTIVI.setup_mudata(mdata, rna_layer="counts", atac_layer="counts")

# PeakVI
scvi.model.PEAKVI.setup_anndata(adata, batch_key="batch")

# veloVI
scvi.external.VELOVI.setup_anndata(adata, spliced_layer="spliced", unspliced_layer="unspliced")
```

### 최소 버전 요구사항

| 패키지 | 최소 버전 | 비고 |
|--------|----------|------|
| scvi-tools | 1.0.0 | 현재 API에 필요 |
| scanpy | 1.9.0 | HVG 선택 개선 |
| anndata | 0.9.0 | 개선된 MuData 지원 |
| torch | 2.0.0 | 성능 개선 |
| mudata | 0.2.0 | MultiVI에 필요 |
| scvelo | 0.2.5 | veloVI에 필요 |

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

# Check if using 1.x API
if hasattr(scvi.model.SCVI, 'setup_anndata'):
    print("Using scvi-tools 1.x API")
else:
    print("WARNING: Using deprecated 0.x API - please upgrade")
```

### 알려진 호환성 문제

| 문제 | 영향 받는 버전 | 해결 방법 |
|------|--------------|----------|
| `setup_anndata` 미발견 | scvi-tools < 1.0 | 1.0+ 으로 업그레이드 |
| MuData 오류 | mudata < 0.2 | `pip install mudata>=0.2.0` |
| CUDA 버전 불일치 | 모든 버전 | 해당 CUDA용 PyTorch 재설치 |
| numpy 2.0 문제 | 2024년 초 빌드 | `pip install numpy<2.0` |

### scvi-tools 업그레이드

```bash
# Upgrade to latest
pip install --upgrade scvi-tools

# Upgrade all dependencies
pip install --upgrade scvi-tools scanpy anndata torch

# If you have issues, clean install
pip uninstall scvi-tools
pip cache purge
pip install scvi-tools
```

## 설치 테스트

```python
# Quick test with sample data
import scvi
import scanpy as sc

# Load test dataset
adata = scvi.data.heart_cell_atlas_subsampled()
print(f"Loaded test data: {adata.shape}")

# Setup and create model (quick test)
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="cell_source")
model = scvi.model.SCVI(adata, n_latent=10)
print("Model created successfully")

# Quick training test (1 epoch)
model.train(max_epochs=1)
print("Training works!")
```
