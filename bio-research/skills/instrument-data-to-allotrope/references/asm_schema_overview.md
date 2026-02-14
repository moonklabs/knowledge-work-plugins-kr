# ASM 스키마 개요

Allotrope Simple Model(ASM)은 시맨틱 일관성을 갖춘 실험실 기기 데이터를 표현하기 위한 JSON 기반 표준입니다.

## 핵심 개념

### 구조
ASM은 계층적 문서 구조를 사용합니다:
- **매니페스트(Manifest)** — 온톨로지 및 스키마 연결
- **데이터(Data)** — 기술(technique)별로 구성된 실제 측정 데이터

### 주요 구성 요소

```json
{
  "$asm.manifest": {
    "vocabulary": ["http://purl.allotrope.org/voc/afo/REC/2023/09/"],
    "contexts": ["http://purl.allotrope.org/json-ld/afo-context-REC-2023-09.jsonld"]
  },
  "<technique>-aggregate-document": {
    "device-system-document": { ... },
    "<technique>-document": [
      {
        "measurement-aggregate-document": {
          "measurement-document": [ ... ]
        }
      }
    ]
  }
}
```

## 필수 메타데이터 문서

### data system document
모든 ASM 출력은 반드시 다음 필드를 포함하는 이 문서를 갖추어야 합니다:
- `ASM file identifier`: 출력 파일 이름
- `data system instance identifier`: 시스템 ID 또는 "N/A"
- `file name`: 소스 입력 파일 이름
- `UNC path`: 소스 파일 경로
- `ASM converter name`: 파서 식별자 (예: "allotropy_beckman_coulter_biomek")
- `ASM converter version`: 버전 문자열
- `software name`: 소스 파일을 생성한 기기 소프트웨어

### device system document
모든 ASM 출력은 반드시 다음 필드를 포함하는 이 문서를 갖추어야 합니다:
- `equipment serial number`: 주 기기 일련번호
- `product manufacturer`: 벤더 이름
- `device document`: 하위 구성 요소 배열 (프로브, 포드 등)
  - `device type`: 표준화된 유형 (예: "liquid handler probe head")
  - `device identifier`: 논리적 이름 (예: "Pod1", 일련번호가 아님)
  - `equipment serial number`: 구성 요소 일련번호
  - `product manufacturer`: 구성 요소 벤더

## 사용 가능한 ASM 기술

공식 ASM 저장소에는 **65개의 기술 스키마**가 포함되어 있습니다:

```
absorbance, automated-reactors, balance, bga, binding-affinity, bulk-density,
cell-counting, cell-culture-analyzer, chromatography, code-reader, conductance,
conductivity, disintegration, dsc, dvs, electronic-lab-notebook,
electronic-spectrometry, electrophoresis, flow-cytometry, fluorescence,
foam-height, foam-qualification, fplc, ftir, gas-chromatography, gc-ms, gloss,
hot-tack, impedance, lc-ms, light-obscuration, liquid-chromatography,
loss-on-drying, luminescence, mass-spectrometry, metabolite-analyzer,
multi-analyte-profiling, nephelometry, nmr, optical-imaging, optical-microscopy,
osmolality, oven-kf, pcr, ph, plate-reader, pressure-monitoring, psd, pumping,
raman, rheometry, sem, solution-analyzer, specific-rotation, spectrophotometry,
stirring, surface-area-analysis, tablet-hardness, temperature-monitoring,
tensile-test, thermogravimetric-analysis, titration, ultraviolet-absorbance,
x-ray-powder-diffraction
```

참조: https://gitlab.com/allotrope-public/asm/-/tree/main/json-schemas/adm

## 기술별 주요 ASM 스키마

아래는 자주 사용되는 기술의 세부 사항입니다:

### 세포 계수(Cell Counting)
스키마: `cell-counting/REC/2024/09/cell-counting.schema.json`

주요 필드:
- `viable-cell-density` (cells/mL)
- `viability` (백분율)
- `total-cell-count`
- `dead-cell-count`
- `cell-diameter-distribution-datum`

### 분광광도법(Spectrophotometry, UV-Vis)
스키마: `spectrophotometry/REC/2024/06/spectrophotometry.schema.json`

주요 필드:
- `absorbance` (무차원)
- `wavelength` (nm)
- `transmittance` (백분율)
- `pathlength` (cm)
- `concentration` (단위 포함)

### 플레이트 리더(Plate Reader)
스키마: `plate-reader/REC/2024/06/plate-reader.schema.json`

주요 필드:
- `absorbance`
- `fluorescence`
- `luminescence`
- `well-location` (A1-H12)
- `plate-identifier`

### qPCR
스키마: `pcr/REC/2024/06/pcr.schema.json`

주요 필드:
- `cycle-threshold-result`
- `amplification-efficiency`
- `melt-curve-datum`
- `target-DNA-description`

### 크로마토그래피(Chromatography)
스키마: `liquid-chromatography/REC/2023/09/liquid-chromatography.schema.json`

주요 필드:
- `retention-time` (분)
- `peak-area`
- `peak-height`
- `peak-width`
- `chromatogram-data-cube`

## 데이터 패턴

### 값 데이텀(Value Datum)
단위가 포함된 단순 값:
```json
{
  "value": 1.5,
  "unit": "mL"
}
```

### 집합 데이텀(Aggregate Datum)
관련 값의 컬렉션:
```json
{
  "measurement-aggregate-document": {
    "measurement-document": [
      { "viable-cell-density": {"value": 2.5e6, "unit": "(cell/mL)"} },
      { "viability": {"value": 95.2, "unit": "%"} }
    ]
  }
}
```

### 데이터 큐브(Data Cube)
다차원 배열 데이터:
```json
{
  "cube-structure": {
    "dimensions": [{"@componentDatatype": "double", "concept": "elapsed time"}],
    "measures": [{"@componentDatatype": "double", "concept": "absorbance"}]
  },
  "data": {
    "dimensions": [[0, 1, 2, 3, 4]],
    "measures": [[0.1, 0.2, 0.3, 0.4, 0.5]]
  }
}
```

## 검증

공식 스키마에 대한 ASM 출력 검증:

```python
import json
import jsonschema
from urllib.request import urlopen

# Load ASM output
with open("output.json") as f:
    asm = json.load(f)

# Get schema URL from manifest
schema_url = asm.get("$asm.manifest", {}).get("$ref")

# Validate (simplified - real validation more complex)
# Note: Full validation requires resolving $ref references
```

## 스키마 저장소

공식 스키마: https://gitlab.com/allotrope-public/asm/-/tree/main/json-schemas/adm

스키마 구조:
```
json-schemas/adm/
├── cell-counting/
│   └── REC/2024/09/
│       └── cell-counting.schema.json
├── spectrophotometry/
│   └── REC/2024/06/
│       └── spectrophotometry.schema.json
├── plate-reader/
│   └── REC/2024/06/
│       └── plate-reader.schema.json
└── ...
```

## 일반적인 문제

### 필드 누락
모든 기기 내보내기가 모든 ASM 필드를 포함하지는 않습니다. 완성도를 보고하십시오:
```python
def report_completeness(asm, expected_fields):
    found = set(extract_all_fields(asm))
    missing = expected_fields - found
    return len(found) / len(expected_fields) * 100
```

### 단위 변형
기기마다 다른 단위 형식을 사용할 수 있습니다. allotropy 라이브러리가 이를 정규화합니다:
- "cells/mL" → "(cell/mL)"
- "%" → "%"
- "nm" → "nm"

### 날짜 형식
ASM은 ISO 8601을 사용합니다: `2024-01-15T10:30:00Z`
