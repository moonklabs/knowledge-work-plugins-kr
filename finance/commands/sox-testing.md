---
description: SOX 샘플 선정, 테스트 워크페이퍼, 통제 평가 생성
argument-hint: "<control area> [period]"
---

# SOX Compliance Testing

> 낯선 플레이스홀더가 보이거나 어떤 도구가 연결되어 있는지 확인하려면 [CONNECTORS.md](../CONNECTORS.md)를 보세요.

**중요**: 이 커맨드는 SOX 워크플로를 보조하지만 감사/법률 자문을 제공하지 않습니다. 산출물은 전문가 검토 후 사용해야 합니다.

## Usage

```bash
/sox <control-area> <period>
```

### Arguments

- `control-area`: `revenue-recognition`, `p2p`, `payroll`, `financial-close`, `treasury`, `fixed-assets`, `inventory`, `itgc`, `entity-level`, `journal-entries` 또는 통제 ID
- `period`: `2024-Q4`, `2024`, `2024-H2` 등

## Workflow

### 1. 테스트 대상 통제 식별

통제 매트릭스를 작성:
- 통제 ID/설명
- 유형(수동/자동/IT의존수동)
- 빈도
- Key 여부
- 위험도
- 재무제표 주장(CEAVOP)

### 2. 샘플 사이즈 산정

빈도/위험도/전기 결과/감사 의존도 기준으로 샘플 수를 산정합니다.

### 3. 샘플 선정

- 랜덤(기본)
- 시스템틱(주기 통제)
- 타겟(고위험 보완)

선정 근거와 표본 목록을 문서화합니다.

### 4. 테스트 워크페이퍼 생성

포함 항목:
- 통제 목적
- 테스트 절차
- 기대 증빙
- 샘플별 통과/실패 결과
- 예외사항/원인/영향
- 결론(Effective, Effective with exceptions, Deficiency, Significant Deficiency, Material Weakness)

### 5. 영역별 템플릿 제공

- Revenue recognition
- Procure to pay
- Financial close
- ITGC

### 6. 결함 평가

Deficiency, Significant Deficiency, Material Weakness 기준으로 분류하고 근거를 적습니다.

### 7. 출력

1. 통제 매트릭스
2. 샘플 선정 결과 및 방법론
3. 테스트 워크페이퍼 템플릿
4. 결과 문서 템플릿
5. 결함 평가 프레임워크
6. 개선/시정 조치 제안
