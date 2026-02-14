---
name: instrument-data-to-allotrope
description: 실험실 기기 출력 파일(PDF, CSV, Excel, TXT)을 Allotrope Simple Model(ASM) JSON 형식 또는 평탄화된 2D CSV로 변환합니다. 과학자가 LIMS 시스템, 데이터 레이크 또는 다운스트림 분석을 위해 기기 데이터를 표준화해야 할 때 이 스킬을 사용합니다. 기기 유형 자동 감지를 지원합니다. 출력에는 전체 ASM JSON, 쉬운 가져오기를 위한 평탄화된 CSV, 데이터 엔지니어를 위한 내보내기 가능한 Python 코드가 포함됩니다. 일반적인 트리거에는 기기 파일 변환, 실험실 데이터 표준화, LIMS/ELN 시스템 업로드를 위한 데이터 준비 또는 프로덕션 파이프라인용 파서 코드 생성이 포함됩니다.
---

# 기기 데이터를 Allotrope로 변환

기기 파일을 LIMS 업로드, 데이터 레이크 또는 데이터 엔지니어링 팀 인수인계를 위한 표준화된 Allotrope Simple Model(ASM) 형식으로 변환합니다.

> **참고: 이것은 예제 스킬입니다**
>
> 이 스킬은 스킬이 데이터 엔지니어링 작업을 어떻게 지원할 수 있는지 보여줍니다—스키마 변환 자동화, 기기 출력 파싱, 프로덕션 준비 코드 생성 등을 수행합니다.
>
> **조직에 맞게 사용자 정의하려면:**
> - `references/` 파일을 수정하여 회사의 특정 스키마 또는 온톨로지 매핑을 포함시키십시오
> - MCP 서버를 사용하여 스키마를 정의하는 시스템(예: LIMS, 데이터 카탈로그 또는 스키마 레지스트리)에 연결하십시오
> - `scripts/`를 확장하여 독점 기기 형식 또는 내부 데이터 표준을 처리하십시오
>
> 이 패턴은 형식 간 변환이나 조직 표준에 대한 검증이 필요한 모든 데이터 변환 워크플로우에 적용할 수 있습니다.

## 워크플로우 개요

1. **기기 유형 감지** — 파일 내용에서 자동 감지 또는 사용자 지정
2. **파일 파싱** — allotropy 라이브러리(네이티브) 또는 유연한 폴백 파서 사용
3. **출력 생성**:
   - ASM JSON (전체 시맨틱 구조)
   - 평탄화된 CSV (2D 테이블 형식)
   - Python 파서 코드 (데이터 엔지니어 인수인계용)
4. **전달** — 요약 및 사용 지침과 함께 파일 제공

> **불확실한 경우:** 필드를 ASM에 어떻게 매핑해야 할지 확실하지 않은 경우(예: 원시 데이터인지 계산된 값인지, 장치 설정인지 환경 조건인지), 사용자에게 확인을 요청하십시오. 가이드는 `references/field_classification_guide.md`를 참조하되, 모호한 경우 추측하지 말고 사용자에게 확인하십시오.

## 빠른 시작

```python
# 필수 패키지 먼저 설치
pip install allotropy pandas openpyxl pdfplumber --break-system-packages

# 핵심 변환 로직
from allotropy.parser_factory import Vendor
from allotropy.to_allotrope import allotrope_from_file

# allotropy를 이용한 변환
asm = allotrope_from_file("instrument_data.csv", Vendor.BECKMAN_VI_CELL_BLU)
```

## 출력 형식 선택

**ASM JSON (기본값)** — 온톨로지 URI가 포함된 전체 시맨틱 구조
- 용도: ASM을 기대하는 LIMS 시스템, 데이터 레이크, 장기 보관
- Allotrope 스키마에 대해 검증

**평탄화된 CSV** — 2D 테이블 표현
- 용도: 빠른 분석, Excel 사용자, JSON을 지원하지 않는 시스템
- 각 측정이 메타데이터가 반복된 하나의 행이 됨

**둘 다** — 최대 유연성을 위해 두 형식 모두 생성

## 계산된 데이터(Calculated Data) 처리

**중요:** 원시 측정값과 계산된/파생된 값을 분리하십시오.

- **원시 데이터** → `measurement-document` (직접 기기 판독값)
- **계산된 데이터** → `calculated-data-aggregate-document` (파생된 값)

계산된 값은 반드시 `data-source-aggregate-document`를 통한 추적성을 포함해야 합니다:

```json
"calculated-data-aggregate-document": {
  "calculated-data-document": [{
    "calculated-data-identifier": "SAMPLE_B1_DIN_001",
    "calculated-data-name": "DNA integrity number",
    "calculated-result": {"value": 9.5, "unit": "(unitless)"},
    "data-source-aggregate-document": {
      "data-source-document": [{
        "data-source-identifier": "SAMPLE_B1_MEASUREMENT",
        "data-source-feature": "electrophoresis trace"
      }]
    }
  }]
}
```

**기기 유형별 일반 계산 필드:**
| 기기 | 계산 필드 |
|------|----------|
| 세포 카운터 | 생존율(Viability) %, 세포 밀도(cell density) 희석 보정값 |
| 분광광도계 | 농도(흡광도에서 산출), 260/280 비율 |
| 플레이트 리더 | 표준 곡선에서의 농도, %CV |
| 전기영동 | DIN/RIN, 영역 농도, 평균 크기 |
| qPCR | 상대 정량값, 배수 변화(fold change) |

원시 데이터와 계산 데이터 분류에 대한 자세한 가이드는 `references/field_classification_guide.md`를 참조하십시오.

## 검증

ASM 출력을 사용자에게 전달하기 전에 항상 검증하십시오:

```bash
python scripts/validate_asm.py output.json
python scripts/validate_asm.py output.json --reference known_good.json  # Compare to reference
python scripts/validate_asm.py output.json --strict  # Treat warnings as errors
```

**검증 규칙:**
- Allotrope ASM 사양 기반 (2024년 12월)
- 최종 업데이트: 2026-01-07
- 소스: https://gitlab.com/allotrope-public/asm

**소프트 검증 방식:**
알 수 없는 기술(technique), 단위 또는 시료 역할은 전방 호환성을 위해 오류가 아닌 **경고**를 생성합니다. Allotrope가 2024년 12월 이후 새 값을 추가하면 검증기가 차단하지 않고 수동 확인을 위해 플래그를 지정합니다. 더 엄격한 검증이 필요한 경우 `--strict` 모드를 사용하여 경고를 오류로 처리하십시오.

**검사 항목:**
- 올바른 기술 선택 (예: 다중 분석물 프로파일링 vs 플레이트 리더)
- 필드 명명 규칙 (하이픈이 아닌 공백으로 구분)
- 계산된 데이터의 추적성 (`data-source-aggregate-document`)
- 측정값 및 계산값에 대한 고유 식별자 존재
- 필수 메타데이터 존재
- 유효한 단위 및 시료 역할 (알 수 없는 값에 대한 소프트 검증)

## 지원 기기

전체 목록은 `references/supported_instruments.md`를 참조하십시오. 주요 기기:

| 카테고리 | 기기 |
|----------|------|
| 세포 계수 | Vi-CELL BLU, Vi-CELL XR, NucleoCounter |
| 분광광도법 | NanoDrop One/Eight/8000, Lunatic |
| 플레이트 리더 | SoftMax Pro, EnVision, Gen5, CLARIOstar |
| ELISA | SoftMax Pro, BMG MARS, MSD Workbench |
| qPCR | QuantStudio, Bio-Rad CFX |
| 크로마토그래피 | Empower, Chromeleon |

## 감지 및 파싱 전략

### 1단계: 네이티브 allotropy 파싱 (권장)
**항상 allotropy를 먼저 시도하십시오.** 사용 가능한 벤더를 직접 확인합니다:

```python
from allotropy.parser_factory import Vendor

# 지원되는 모든 벤더 목록 출력
for v in Vendor:
    print(f"{v.name}")

# 일반적인 벤더 예시:
# AGILENT_TAPESTATION_ANALYSIS  (TapeStation XML용)
# BECKMAN_VI_CELL_BLU
# THERMO_FISHER_NANODROP_EIGHT
# MOLDEV_SOFTMAX_PRO
# APPBIO_QUANTSTUDIO
# ... 기타 다수
```

**사용자가 파일을 제공하면, 수동 파싱으로 폴백하기 전에 allotropy가 지원하는지 확인하십시오.** `scripts/convert_to_asm.py`의 자동 감지는 allotropy 벤더의 일부만 포함합니다.

### 2단계: 유연한 폴백 파싱
**allotropy가 기기를 지원하지 않는 경우에만 사용하십시오.** 이 폴백은:
- `calculated-data-aggregate-document`를 생성하지 않음
- 전체 추적성을 포함하지 않음
- 단순화된 ASM 구조 생성

유연한 파서 사용 기능:
- 열 이름 퍼지 매칭
- 헤더에서 단위 추출
- 파일 구조에서 메타데이터 추출

### 3단계: PDF 추출
PDF 전용 파일의 경우 pdfplumber를 사용하여 테이블을 추출한 후 2단계 파싱을 적용합니다.

## 파싱 전 체크리스트

사용자 정의 파서를 작성하기 전에 반드시 다음을 수행하십시오:

1. **allotropy가 지원하는지 확인** — 사용 가능한 경우 네이티브 파서 사용
2. **참조 ASM 파일 찾기** — `references/examples/`를 확인하거나 사용자에게 요청
3. **기기별 가이드 확인** — `references/instrument_guides/` 확인
4. **참조 파일에 대해 검증** — `validate_asm.py --reference <file>` 실행

## 흔한 실수 방지

| 실수 | 올바른 접근 방식 |
|------|----------------|
| 매니페스트를 객체로 사용 | URL 문자열 사용 |
| 소문자 감지 유형 | "absorbance"가 아닌 "Absorbance" 사용 |
| "emission wavelength setting" | 방출에는 "detector wavelength setting" 사용 |
| 모든 측정을 하나의 문서에 포함 | 웰/시료 위치별로 그룹화 |
| 프로시저 메타데이터 누락 | 측정당 모든 장치 설정 추출 |

## 데이터 엔지니어를 위한 코드 내보내기

과학자가 인수인계할 수 있는 독립형 Python 스크립트를 생성합니다:

```python
# Export parser code
python scripts/export_parser.py --input "data.csv" --vendor "VI_CELL_BLU" --output "parser_script.py"
```

내보낸 스크립트 특징:
- pandas/allotropy 외에 외부 의존성 없음
- 인라인 문서 포함
- Jupyter 노트북에서 실행 가능
- 데이터 파이프라인용 프로덕션 준비 완료

## 파일 구조

```
instrument-data-to-allotrope/
├── SKILL.md                          # This file
├── scripts/
│   ├── convert_to_asm.py            # Main conversion script
│   ├── flatten_asm.py               # ASM → 2D CSV conversion
│   ├── export_parser.py             # Generate standalone parser code
│   └── validate_asm.py              # Validate ASM output quality
└── references/
    ├── supported_instruments.md     # Full instrument list with Vendor enums
    ├── asm_schema_overview.md       # ASM structure reference
    ├── field_classification_guide.md # Where to put different field types
    └── flattening_guide.md          # How flattening works
```

## 사용 예제

### 예제 1: Vi-CELL BLU 파일
```
사용자: "이 세포 계수 데이터를 Allotrope 형식으로 변환해 줘."
[viCell_Results.xlsx 업로드]

Claude:
1. Vi-CELL BLU 감지 (95% 신뢰도)
2. allotropy 네이티브 파서를 사용하여 변환
3. 출력 파일:
   - viCell_Results_asm.json (전체 ASM JSON)
   - viCell_Results_flat.csv (2D 평탄화 형식)
   - viCell_parser.py (내보내기 가능한 Python 코드)
```

### 예제 2: 코드 인수인계 요청
```
사용자: "우리 데이터 엔지니어에게 NanoDrop 파일을 파싱할 수 있는 코드를 전달해야 해."

Claude:
1. 독립적으로 실행 가능한 Python 스크립트 생성
2. 샘플 입력/출력 포함
3. 모든 가정 사항 문서화
4. Jupyter 노트북 버전 제공
```

### 예제 3: LIMS 준비 평탄화 출력
```
사용자: "이 ELISA 데이터를 우리 LIMS에 업로드할 수 있는 CSV로 변환해 줘."

Claude:
1. 플레이트 리더 데이터 파싱
2. 다음 컬럼을 포함하는 평탄화된 CSV 생성:
   - sample_identifier, well_position, measurement_value, measurement_unit
   - instrument_serial_number, analysis_datetime, assay_type
3. 일반적인 LIMS 가져오기 요구 사항에 맞게 검증
```

## 구현 참고사항

### allotropy 설치
```bash
pip install allotropy --break-system-packages
```

### 파싱 실패 처리
allotropy 네이티브 파싱이 실패하는 경우:
1. 디버깅을 위해 오류 로깅
2. 유연한 파서로 폴백
3. 사용자에게 메타데이터 완성도 감소 보고
4. 기기에서 다른 형식으로 내보내기를 제안

### ASM 스키마 검증
사용 가능한 경우 Allotrope 스키마에 대해 출력을 검증합니다:
```python
import jsonschema
# Schema URLs in references/asm_schema_overview.md
```
