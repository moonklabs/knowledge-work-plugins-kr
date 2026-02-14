---
name: single-cell-rna-qc
description: scverse 모범 사례에 따라 단일세포 RNA 시퀀싱 데이터(.h5ad 또는 .h5 파일)에 대한 품질 관리를 수행하며, MAD 기반 필터링과 종합적인 시각화를 제공합니다. 사용자가 QC 분석, 저품질 세포 필터링, 데이터 품질 평가 또는 scverse/scanpy 모범 사례를 요청할 때 사용합니다.
---

# 단일세포 RNA 시퀀싱 품질 관리

scverse 모범 사례에 따른 단일세포 RNA 시퀀싱 데이터의 자동화된 QC 워크플로우입니다.

## 스킬 사용 시점

다음과 같은 요청 시 사용합니다:
- 단일세포 RNA 시퀀싱 데이터에 대한 품질 관리 또는 QC 요청
- 저품질 세포 필터링 또는 데이터 품질 평가
- QC 시각화 또는 지표 필요
- scverse/scanpy 모범 사례 준수 요청
- MAD 기반 필터링 또는 이상치 탐지 요청

**지원 입력 형식:**
- `.h5ad` 파일 (scanpy/Python 워크플로우의 AnnData 형식)
- `.h5` 파일 (10X Genomics Cell Ranger 출력)

**기본 권장 사항**: 사용자에게 특별한 맞춤 요구 사항이 있거나 비표준 필터링 로직을 명시적으로 요청하지 않는 한 접근법 1(전체 파이프라인)을 사용하십시오.

## 접근법 1: 전체 QC 파이프라인 (표준 워크플로우에 권장)

scverse 모범 사례에 따른 표준 QC의 경우, 편의 스크립트 `scripts/qc_analysis.py`를 사용합니다:

```bash
python3 scripts/qc_analysis.py input.h5ad
# or for 10X Genomics .h5 files:
python3 scripts/qc_analysis.py raw_feature_bc_matrix.h5
```

이 스크립트는 파일 형식을 자동으로 감지하여 적절하게 로드합니다.

**사용 시점:**
- 조정 가능한 임계값을 사용한 표준 QC 워크플로우 (모든 세포에 동일한 필터링 적용)
- 다수의 데이터셋 일괄 처리
- 빠른 탐색적 분석
- "바로 작동하는" 솔루션이 필요한 경우

**요구 사항:** anndata, scanpy, scipy, matplotlib, seaborn, numpy

**파라미터:**

명령줄 파라미터를 사용하여 필터링 임계값과 유전자 패턴을 맞춤 설정할 수 있습니다:
- `--output-dir` - 출력 디렉토리
- `--mad-counts`, `--mad-genes`, `--mad-mt` - 카운트/유전자/MT%에 대한 MAD 임계값
- `--mt-threshold` - 미토콘드리아 % 고정 임계값
- `--min-cells` - 유전자 필터링 임계값
- `--mt-pattern`, `--ribo-pattern`, `--hb-pattern` - 종별 유전자 이름 패턴

`--help`를 사용하면 현재 기본값을 확인할 수 있습니다.

**출력:**

모든 파일은 기본적으로 `<input_basename>_qc_results/` 디렉토리에 저장됩니다(`--output-dir`로 지정된 디렉토리에 저장 가능):
- `qc_metrics_before_filtering.png` - 필터링 전 시각화
- `qc_filtering_thresholds.png` - MAD 기반 임계값 오버레이
- `qc_metrics_after_filtering.png` - 필터링 후 품질 지표
- `<input_basename>_filtered.h5ad` - 다운스트림 분석에 사용 가능한 정제된 필터링 데이터셋
- `<input_basename>_with_qc.h5ad` - QC 주석이 보존된 원본 데이터

사용자 접근을 위해 출력을 복사할 때는 전체 디렉토리가 아닌 개별 파일을 복사하여 사용자가 직접 미리볼 수 있도록 합니다.

### 워크플로우 단계

이 스크립트는 다음 단계를 수행합니다:

1. **QC 지표 계산** - 카운트 깊이, 유전자 탐지, 미토콘드리아/리보솜/헤모글로빈 함량
2. **MAD 기반 필터링 적용** - MAD 임계값을 사용한 허용적 이상치 탐지 (카운트/유전자/MT%)
3. **유전자 필터링** - 소수의 세포에서만 탐지된 유전자 제거
4. **시각화 생성** - 임계값 오버레이를 포함한 종합적인 전후 플롯

## 접근법 2: 모듈식 빌딩 블록 (맞춤 워크플로우용)

맞춤 분석 워크플로우 또는 비표준 요구 사항의 경우, `scripts/qc_core.py` 및 `scripts/qc_plotting.py`의 모듈식 유틸리티 함수를 사용합니다:

```python
# Run from scripts/ directory, or add scripts/ to sys.path if needed
import anndata as ad
from qc_core import calculate_qc_metrics, detect_outliers_mad, filter_cells
from qc_plotting import plot_qc_distributions  # Only if visualization needed

adata = ad.read_h5ad('input.h5ad')
calculate_qc_metrics(adata, inplace=True)
# ... custom analysis logic here
```

**사용 시점:**
- 다른 워크플로우 필요 시 (단계 건너뛰기, 순서 변경, 하위 집합에 다른 임계값 적용)
- 조건부 로직 (예: 뉴런과 다른 세포에 다른 필터링 적용)
- 부분 실행 (지표/시각화만 수행, 필터링 없음)
- 더 큰 파이프라인의 다른 분석 단계와 통합
- 명령줄 파라미터로 지원되지 않는 맞춤 필터링 기준

**사용 가능한 유틸리티 함수:**

`qc_core.py` (핵심 QC 연산):
- `calculate_qc_metrics(adata, mt_pattern, ribo_pattern, hb_pattern, inplace=True)` - QC 지표 계산 및 adata 주석
- `detect_outliers_mad(adata, metric, n_mads, verbose=True)` - MAD 기반 이상치 탐지, 불리언 마스크 반환
- `apply_hard_threshold(adata, metric, threshold, operator='>', verbose=True)` - 고정 임계값 적용, 불리언 마스크 반환
- `filter_cells(adata, mask, inplace=False)` - 불리언 마스크를 적용하여 세포 필터링
- `filter_genes(adata, min_cells=20, min_counts=None, inplace=True)` - 탐지 기준으로 유전자 필터링
- `print_qc_summary(adata, label='')` - 요약 통계 출력

`qc_plotting.py` (시각화):
- `plot_qc_distributions(adata, output_path, title)` - 종합 QC 플롯 생성
- `plot_filtering_thresholds(adata, outlier_masks, thresholds, output_path)` - 필터링 임계값 시각화
- `plot_qc_after_filtering(adata, output_path)` - 필터링 후 플롯 생성

**맞춤 워크플로우 예시:**

**예시 1: 지표 계산 및 시각화만 수행, 필터링은 아직 하지 않음**
```python
adata = ad.read_h5ad('input.h5ad')
calculate_qc_metrics(adata, inplace=True)
plot_qc_distributions(adata, 'qc_before.png', title='Initial QC')
print_qc_summary(adata, label='Before filtering')
```

**예시 2: MT% 필터링만 적용, 다른 지표는 허용적으로 유지**
```python
adata = ad.read_h5ad('input.h5ad')
calculate_qc_metrics(adata, inplace=True)

# Only filter high MT% cells
high_mt = apply_hard_threshold(adata, 'pct_counts_mt', 10, operator='>')
adata_filtered = filter_cells(adata, ~high_mt)
adata_filtered.write('filtered.h5ad')
```

**예시 3: 하위 집합별 다른 임계값 적용**
```python
adata = ad.read_h5ad('input.h5ad')
calculate_qc_metrics(adata, inplace=True)

# Apply type-specific QC (assumes cell_type metadata exists)
neurons = adata.obs['cell_type'] == 'neuron'
other_cells = ~neurons

# Neurons tolerate higher MT%, other cells use stricter threshold
neuron_qc = apply_hard_threshold(adata[neurons], 'pct_counts_mt', 15, operator='>')
other_qc = apply_hard_threshold(adata[other_cells], 'pct_counts_mt', 8, operator='>')
```

## 모범 사례

1. **필터링은 허용적으로** - 기본 임계값은 희귀 세포 집단의 손실을 방지하기 위해 의도적으로 대부분의 세포를 유지합니다
2. **시각화 검토** - 필터링이 생물학적으로 타당한지 확인하기 위해 항상 전후 플롯을 검토하십시오
3. **데이터셋별 요인 고려** - 일부 조직은 자연적으로 미토콘드리아 함량이 높습니다 (예: 뉴런, 심근세포)
4. **유전자 주석 확인** - 미토콘드리아 유전자 접두사는 종에 따라 다릅니다 (마우스: mt-, 인간: MT-)
5. **필요시 반복** - QC 파라미터는 특정 실험 또는 조직 유형에 따라 조정이 필요할 수 있습니다

## 참고 자료

자세한 QC 방법론, 파라미터 근거 및 문제 해결 지침은 `references/scverse_qc_guidelines.md`를 참조하십시오. 이 레퍼런스는 다음을 제공합니다:
- 각 QC 지표의 상세 설명 및 중요성
- MAD 기반 임계값의 근거 및 고정 임계값보다 나은 이유
- QC 시각화 해석 지침 (히스토그램, 바이올린 플롯, 산점도)
- 종별 유전자 주석 고려 사항
- 필터링 파라미터 조정 시기 및 방법
- 고급 QC 고려 사항 (주변 RNA 보정, 더블렛 탐지)

방법론에 대한 심층적인 이해가 필요하거나 QC 문제를 해결할 때 이 레퍼런스를 로드하십시오.

## QC 이후 다음 단계

일반적인 다운스트림 분석 단계:
- 주변 RNA 보정 (SoupX, CellBender)
- 더블렛 탐지 (scDblFinder)
- 정규화 (log-normalize, scran)
- 피처 선택 및 차원 축소
- 클러스터링 및 세포 유형 주석
