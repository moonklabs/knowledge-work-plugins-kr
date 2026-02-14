# 필드 분류 가이드

이 가이드는 기기 데이터 필드를 올바른 ASM 문서 위치로 분류하는 데 도움을 줍니다. 원시 기기 출력을 Allotrope Simple Model 구조에 매핑할 때 사용하십시오.

## ASM 문서 계층 구조

```
<technique>-aggregate-document
├── device-system-document          # 기기 하드웨어 정보
├── data-system-document            # 소프트웨어/변환 정보
├── <technique>-document[]          # 실행/시퀀스별 데이터
│   ├── analyst                     # 분석 수행자
│   ├── measurement-aggregate-document
│   │   ├── measurement-time
│   │   ├── measurement-document[]  # 개별 측정값
│   │   │   ├── sample-document
│   │   │   ├── device-control-aggregate-document
│   │   │   └── [measurement fields]
│   │   └── [aggregate-level metadata]
│   ├── processed-data-aggregate-document
│   │   └── processed-data-document[]
│   │       ├── data-processing-document
│   │       └── [processed results]
│   └── calculated-data-aggregate-document
│       └── calculated-data-document[]
```

## 필드 분류 카테고리

### 1. 장치/기기 정보 → `device-system-document`

물리적 기기의 하드웨어 및 펌웨어 세부 정보입니다.

| 필드 유형 | ASM 필드 | 예시 |
|-----------|----------|------|
| 기기 이름 | `model-number` | "Vi-CELL BLU", "NanoDrop One" |
| 일련번호 | `equipment-serial-number` | "VCB-12345", "SN001234" |
| 제조사 | `product-manufacturer` | "Beckman Coulter", "Thermo Fisher" |
| 펌웨어 버전 | `firmware-version` | "v2.1.3" |
| 장치 ID | `device-identifier` | "Instrument_01" |
| 브랜드 | `brand-name` | "Beckman Coulter" |

**규칙:** 값이 물리적 기기를 설명하며 실행 간에 변경되지 않으면 `device-system-document`에 배치합니다.

---

### 2. 소프트웨어/데이터 시스템 정보 → `data-system-document`

데이터 수집, 분석 또는 변환에 사용된 소프트웨어 정보입니다.

| 필드 유형 | ASM 필드 | 예시 |
|-----------|----------|------|
| 소프트웨어 이름 | `software-name` | "Chromeleon", "Gen5" |
| 소프트웨어 버전 | `software-version` | "7.3.2" |
| 파일 이름 | `file-name` | "experiment_001.xlsx" |
| 파일 경로 | `file-identifier` | "/data/runs/2024-01-15/" |
| 데이터베이스 ID | `ASM-converter-name` | "allotropy v0.1.55" |

**규칙:** 값이 소프트웨어, 파일 메타데이터 또는 데이터 출처를 설명하면 `data-system-document`에 배치합니다.

---

### 3. 시료 정보 → `sample-document`

분석 대상 생물학적/화학적 시료에 대한 메타데이터입니다.

| 필드 유형 | ASM 필드 | 예시 |
|-----------|----------|------|
| 시료 ID | `sample-identifier` | "Sample_A", "LIMS-001234" |
| 시료 이름 | `written-name` | "CHO Cell Culture Day 5" |
| 시료 유형/역할 | `sample-role-type` | "unknown sample role", "control sample role" |
| 배치 ID | `batch-identifier` | "Batch-2024-001" |
| 설명 | `description` | "Protein expression sample" |
| 웰 위치 | `location-identifier` | "A1", "B3" |

**규칙:** 값이 무엇을 측정했는지(어떻게가 아닌)를 식별하거나 설명하면 `sample-document`에 배치합니다.

---

### 4. 장치 제어 설정 → `device-control-aggregate-document`

측정 중 사용된 기기 설정 및 매개변수입니다.

| 필드 유형 | ASM 필드 | 예시 |
|-----------|----------|------|
| 주입 부피 | `sample-volume-setting` | 10 µL |
| 파장 | `detector-wavelength-setting` | 254 nm |
| 온도 | `compartment-temperature` | 37°C |
| 유속 | `flow-rate` | 1.0 mL/min |
| 노출 시간 | `exposure-duration-setting` | 500 ms |
| 검출기 게인 | `detector-gain-setting` | 1.5 |
| 조명 | `illumination-setting` | 80% |

**규칙:** 값이 측정에 영향을 미치는 구성 가능한 기기 매개변수이면 `device-control-aggregate-document`에 배치합니다.

---

### 5. 환경 조건 → `device-control-document` 또는 기술별 위치

측정 중 주변 또는 제어된 환경 매개변수입니다.

| 필드 유형 | ASM 필드 | 예시 |
|-----------|----------|------|
| 주변 온도 | `ambient-temperature` | 22.5°C |
| 습도 | `ambient-relative-humidity` | 45% |
| 컬럼 온도 | `compartment-temperature` | 30°C |
| 시료 온도 | `sample-temperature` | 4°C |
| 전기영동 온도 | (기술별) | 26.4°C |

**규칙:** 측정 품질에 영향을 미치는 환경 조건은 장치 제어 또는 기술별 위치에 배치합니다.

---

### 6. 원시 측정 데이터 → `measurement-document`

직접 기기 판독값 — "근거 데이터(ground truth)"입니다.

| 필드 유형 | ASM 필드 | 예시 |
|-----------|----------|------|
| 흡광도 | `absorbance` | 0.523 AU |
| 형광 | `fluorescence` | 12500 RFU |
| 세포 수 | `total-cell-count` | 2.5e6 cells |
| 피크 면적 | `peak-area` | 1234.5 mAU·min |
| 머무름 시간 | `retention-time` | 5.67 min |
| Ct 값 | `cycle-threshold-result` | 24.5 |
| 농도 (측정값) | `mass-concentration` | 1.5 mg/mL |

**규칙:** 값이 이 분석에서 다른 값으로부터 계산되지 않은 직접 기기 판독값이면 `measurement-document`에 배치합니다.

---

### 7. 계산/파생 데이터 → `calculated-data-aggregate-document`

원시 측정값에서 계산된 값입니다.

| 필드 유형 | ASM 필드 | 예시 |
|-----------|----------|------|
| 생존율 % | `calculated-result` | 95.2% |
| 농도 (표준 곡선에서) | `calculated-result` | 125 ng/µL |
| 비율 (260/280) | `calculated-result` | 1.89 |
| 상대 정량값 | `calculated-result` | 2.5x |
| % 회수율 | `calculated-result` | 98.7% |
| CV% | `calculated-result` | 2.3% |

**계산 데이터 문서 구조:**
```json
{
  "calculated-data-name": "viability",
  "calculated-result": {"value": 95.2, "unit": "%"},
  "calculation-description": "viable cells / total cells * 100"
}
```

**규칙:** 값이 이 분석의 다른 측정값에서 계산된 것이면 `calculated-data-aggregate-document`에 배치합니다. 가능하면 `calculation-description`을 포함하십시오.

---

### 8. 처리/분석 데이터 → `processed-data-aggregate-document`

데이터 처리 알고리즘(피크 적분, 세포 분류 등)의 결과입니다.

| 필드 유형 | ASM 필드 | 예시 |
|-----------|----------|------|
| 피크 목록 | `peak-list` | 적분된 피크 결과 |
| 세포 크기 분포 | `cell-diameter-distribution` | 히스토그램 데이터 |
| 기선 보정 데이터 | (processed-data-document 내) | 보정된 스펙트럼 |
| 피팅 곡선 | (processed-data-document 내) | 표준 곡선 피팅 |

**관련 `data-processing-document`:**
```json
{
  "cell-type-processing-method": "trypan blue exclusion",
  "cell-density-dilution-factor": {"value": 2, "unit": "(unitless)"},
  "minimum-cell-diameter-setting": {"value": 5, "unit": "µm"},
  "maximum-cell-diameter-setting": {"value": 50, "unit": "µm"}
}
```

**규칙:** 값이 원시 데이터에 적용된 알고리즘 또는 처리 방법의 결과이면 `processed-data-aggregate-document`에 배치하고, 처리 매개변수는 `data-processing-document`에 포함합니다.

---

### 9. 타이밍/타임스탬프 → 다양한 위치

| 타임스탬프 유형 | 위치 | ASM 필드 |
|----------------|------|----------|
| 측정 시간 | `measurement-document` | `measurement-time` |
| 실행 시작 시간 | `analysis-sequence-document` | `analysis-sequence-start-time` |
| 실행 종료 시간 | `analysis-sequence-document` | `analysis-sequence-end-time` |
| 데이터 내보내기 시간 | `data-system-document` | (사용자 정의) |

**규칙:** ISO 8601 형식을 사용합니다: `2024-01-15T10:30:00Z`

---

### 10. 분석자/작업자 정보 → `<technique>-document`

| 필드 유형 | ASM 필드 | 예시 |
|-----------|----------|------|
| 작업자 이름 | `analyst` | "jsmith" |
| 검토자 | (사용자 정의 또는 확장) | "Pending" |

**규칙:** 분석자는 개별 측정이 아닌 technique-document 수준에 배치합니다.

---

## 의사결정 트리

```
이 필드가 설명하는 것은...

기기 자체?
├── 하드웨어 사양 → device-system-document
└── 소프트웨어/파일 → data-system-document

시료?
└── 시료 ID, 이름, 유형, 배치 → sample-document

기기 설정?
└── 구성 가능한 매개변수 → device-control-aggregate-document

환경 조건?
└── 온도, 습도 등 → device-control-document

직접 판독값?
└── 원시 기기 출력 → measurement-document

계산된 값?
├── 다른 측정값에서 산출 → calculated-data-document
└── 처리 알고리즘에서 산출 → processed-data-document

시간?
├── 측정 시점 → measurement-document.measurement-time
└── 실행 시작/종료 시점 → analysis-sequence-document

누가 수행했는가?
└── 작업자/분석자 → <technique>-document.analyst
```

## 일반적인 기기-ASM 매핑

> **참고:** 이 매핑은 [Benchling allotropy 라이브러리](https://github.com/Benchling-Open-Source/allotropy/tree/main/src/allotropy/parsers)에서 파생되었습니다. 정확한 매핑은 특정 기기의 파서 소스 코드를 참조하십시오.

### 세포 카운터 (Vi-CELL BLU)
*소스: `allotropy/parsers/beckman_vi_cell_blu/vi_cell_blu_structure.py`*

| 기기 필드 | ASM 필드 |
|-----------|----------|
| Sample ID | `sample_identifier` |
| Analysis date/time | `measurement_time` |
| Analysis by | `analyst` |
| Viability (%) | `viability` |
| Viable (x10^6) cells/mL | `viable_cell_density` |
| Total (x10^6) cells/mL | `total_cell_density` |
| Cell count | `total_cell_count` |
| Viable cells | `viable_cell_count` |
| Average diameter (μm) | `average_total_cell_diameter` |
| Average viable diameter (μm) | `average_live_cell_diameter` |
| Average circularity | `average_total_cell_circularity` |
| Cell type | `cell_type_processing_method` (data-processing) |
| Dilution | `cell_density_dilution_factor` (data-processing) |
| Min/Max Diameter | `minimum/maximum_cell_diameter_setting` (data-processing) |

### 분광광도계 (NanoDrop)
| 기기 필드 | ASM 필드 |
|-----------|----------|
| Sample Name | `sample_identifier` |
| A260, A280 | `absorbance` (파장 포함) |
| Concentration | `mass_concentration` |
| 260/280 ratio | `a260_a280_ratio` |
| Pathlength | `pathlength` |

### 플레이트 리더
| 기기 필드 | ASM 필드 |
|-----------|----------|
| Well | `location_identifier` |
| Sample Type | `sample_role_type` |
| Absorbance/OD | `absorbance` |
| Fluorescence | `fluorescence` |
| Plate ID | `container_identifier` |

### 크로마토그래피 (HPLC)
| 기기 필드 | ASM 필드 |
|-----------|----------|
| Sample ID | `sample_identifier` |
| Injection Volume | `injection_volume` |
| Retention Time | `retention_time` |
| Peak Area | `peak_area` |
| Peak Height | `peak_height` |
| Column Temp | `column_oven_temperature` |
| Flow Rate | `flow_rate` |

## 단위 처리

소스 데이터에 명시적으로 존재하는 단위만 사용하십시오. 값에 단위가 지정되지 않은 경우:
- 단위 값으로 `(unitless)` 사용
- 도메인 지식에 기반하여 단위를 추론하지 마십시오

## 계산 데이터 추적성

계산 값을 생성할 때는 항상 `data-source-aggregate-document`를 사용하여 소스 데이터에 연결하십시오:

```json
{
    "calculated-data-name": "DIN",
    "calculated-result": {"value": 5.8, "unit": "(unitless)"},
    "calculated-data-identifier": "TEST_ID_147",
    "data-source-aggregate-document": {
        "data-source-document": [{
            "data-source-identifier": "TEST_ID_145",
            "data-source-feature": "sample"
        }]
    }
}
```

이는 "DIN 5.8이 `TEST_ID_145`의 시료에서 계산되었다"를 선언합니다.

**이것이 중요한 이유:**
- **감사**: 값이 특정 원시 데이터에서 나왔음을 증명
- **디버깅**: 예상치 못한 결과를 소스까지 추적
- **재처리**: 알고리즘이 변경되면 어떤 입력을 다시 분석할지 파악

**고유 ID를 부여해야 하는 항목:**
- 측정값, 피크, 영역 및 계산 값
- 일관된 명명 패턴 사용 (예: `INSTRUMENT_TYPE_TEST_ID_N`)

이를 통해 양방향 순회가 가능합니다: 계산 → 원시, 또는 원시 → 모든 파생 값을 추적할 수 있습니다.

---

## 중첩 문서 구조 (핵심)

흔한 실수는 중첩 구조로 래핑해야 하는 필드를 측정 문서에 직접 "평탄화"하는 것입니다. 이는 스키마 준수를 깨뜨리고 시맨틱 컨텍스트를 잃습니다.

### 중첩이 중요한 이유

ASM은 시맨틱 그룹화를 위해 중첩 문서를 사용합니다:

| 문서 | 목적 | 포함 내용 |
|------|------|----------|
| `sample document` | 무엇을 측정했는지 | 시료 ID, 위치, 플레이트 식별자 |
| `device control aggregate document` | 기기 작동 방식 | 설정, 매개변수, 기술 |
| `custom information document` | 벤더 특정 필드 | 표준 ASM 용어에 매핑되지 않는 필드 |

### Sample Document 필드

다음 필드는 반드시 `sample document` 내에 있어야 하며, 측정에 평탄화하면 안 됩니다:

```json
// ❌ 잘못됨 - 필드가 측정에 평탄화됨
{
  "measurement identifier": "TEST_001",
  "sample identifier": "Sample_A",
  "location identifier": "A1",
  "absorbance": {"value": 0.5, "unit": "(unitless)"}
}

// ✅ 올바름 - 필드가 sample document에 중첩됨
{
  "measurement identifier": "TEST_001",
  "sample document": {
    "sample identifier": "Sample_A",
    "location identifier": "A1",
    "well plate identifier": "96WP001"
  },
  "absorbance": {"value": 0.5, "unit": "(unitless)"}
}
```

**sample document에 속하는 필드:**
- `sample identifier` — 시료 ID/이름
- `written name` — 설명적 시료 이름
- `batch identifier` — 배치/로트 번호
- `sample role type` — 표준, 블랭크, 대조군, 미지 시료
- `location identifier` — 웰 위치 (A1, B3 등)
- `well plate identifier` — 플레이트 바코드
- `description` — 시료 설명

### Device Control Document 필드

기기 설정은 반드시 `device control aggregate document` 내에 있어야 합니다:

```json
// ❌ 잘못됨 - 장치 설정이 평탄화됨
{
  "measurement identifier": "TEST_001",
  "device identifier": "Pod1",
  "technique": "Custom",
  "volume": {"value": 26, "unit": "μL"}
}

// ✅ 올바름 - 설정이 device control에 중첩됨
{
  "measurement identifier": "TEST_001",
  "device control aggregate document": {
    "device control document": [{
      "device type": "liquid handler",
      "device identifier": "Pod1"
    }]
  },
  "aspiration volume": {"value": 26, "unit": "μL"}
}
```

**device control에 속하는 필드:**
- `device type` — 장치 유형
- `device identifier` — 장치 ID
- `detector wavelength setting` — 검출 파장
- `compartment temperature` — 온도 설정
- `sample volume setting` — 부피 설정
- `flow rate` — 유속 설정

### Custom Information Document

표준 ASM 용어에 매핑되지 않는 벤더 특정 필드는 `custom information document`에 배치합니다:

```json
"device control document": [{
  "device type": "liquid handler",
  "custom information document": {
    "probe": "2",
    "pod": "Pod1",
    "source labware name": "Inducer",
    "destination labware name": "GRP1"
  }
}]
```

### 액체 핸들러: 이송 페어링

액체 핸들러의 경우 측정은 별도 작업이 아닌 완전한 이송(흡인 + 분주)을 나타냅니다:

```json
// ❌ 잘못됨 - 흡인과 분주를 별도 레코드로 분리
[
  {"measurement identifier": "OP_001", "transfer type": "Aspirate", "volume": {"value": 26, "unit": "μL"}},
  {"measurement identifier": "OP_002", "transfer type": "Dispense", "volume": {"value": 26, "unit": "μL"}}
]

// ✅ 올바름 - 소스와 대상이 포함된 단일 레코드
{
  "measurement identifier": "TRANSFER_001",
  "sample document": {
    "source well location identifier": "1",
    "destination well location identifier": "2",
    "source well plate identifier": "96WP001",
    "destination well plate identifier": "96WP002"
  },
  "aspiration volume": {"value": 26, "unit": "μL"},
  "transfer volume": {"value": 26, "unit": "μL"}
}
```

**페어링 로직:**
1. 프로브 번호로 흡인과 분주 작업 매칭
2. 매칭된 쌍당 하나의 측정 생성
3. 흡인 위치에는 `source_*` 필드 사용
4. 분주 위치에는 `destination_*` 필드 사용
5. `aspiration volume`과 `transfer volume` 모두 포함

### 빠른 참조: 중첩 결정

```
이 필드가 설명하는 것은...

측정 대상 시료?
├── 시료 ID, 이름, 배치 → sample document
├── 웰 위치 → sample document.location identifier
├── 플레이트 바코드 → sample document.well plate identifier
└── 소스/대상 위치 → sample document (접두사 포함)

기기 설정?
├── 표준 설정 → device control aggregate document
└── 벤더 특정 → custom information document

측정값?
└── 측정 문서에 직접 배치 (예: absorbance, volume)

이송 작업 유형?
└── "transfer type"을 사용하지 마십시오 - 대신 소스/대상 필드가 있는
    단일 측정으로 페어링
```

### 검증

중첩 문제를 확인하려면 `validate_asm.py`를 사용하십시오:
```bash
python scripts/validate_asm.py output.json --reference known_good.json
```

검증기가 확인하는 항목:
- 측정에 잘못 평탄화된 필드
- 누락된 `sample document` 래퍼
- 누락된 `device control aggregate document` 래퍼
- 벤더 필드에 대한 누락된 `custom information document`
- 액체 핸들러: 페어링된 레코드 대신 별도의 이송 유형

## 참고 자료

- [Allotrope Simple Model 소개](https://www.allotrope.org/introduction-to-allotrope-simple-model)
- [Benchling allotropy 라이브러리](https://github.com/Benchling-Open-Source/allotropy)
- [Allotrope Foundation ASM 개요](https://www.allotrope.org/asm)
