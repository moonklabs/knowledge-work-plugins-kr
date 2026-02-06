---
description: 조직 협상 플레이북 기준 계약서 검토, 이탈 항목 플래그, 레드라인 제안, 비즈니스 영향 분석
argument-hint: "<contract file or text>"
---

# /review-contract -- Contract Review Against Playbook

> 낯선 플레이스홀더가 보이거나 어떤 도구가 연결되어 있는지 확인하려면 [CONNECTORS.md](../CONNECTORS.md)를 보세요.

조직의 협상 플레이북 기준으로 계약서를 조항별 검토하고 이탈 항목, 레드라인, 비즈니스 영향을 제공합니다.

## Invocation

```bash
/review-contract
```

## Workflow

### 1. 계약서 입력

입력 형식:
- 파일 업로드(PDF, DOCX 등)
- URL(CLM, 문서 저장소)
- 텍스트 붙여넣기

### 2. 컨텍스트 수집

- 우리 측 지위(벤더/고객 등)
- 마감 일정
- 중점 검토 영역
- 딜 맥락(규모, 전략성, 기존 관계)

### 3. 플레이북 로드

`legal.local.md` 등에서 표준 입장/허용 범위/에스컬레이션 트리거를 로드합니다.
없으면:
1. 플레이북 설정 지원
2. 일반 상관행 기준으로 임시 검토(한계 명시)

### 4. 조항별 분석

최소 검토 범주:
- 책임 제한
- 면책
- IP 소유권
- 데이터 보호
- 비밀유지
- 진술/보증
- 계약기간/종료
- 준거법/분쟁해결
- 보험
- 양도
- 불가항력
- 지급조건

### 5. 이탈 분류

- `GREEN`: 수용 가능
- `YELLOW`: 협상 필요
- `RED`: 에스컬레이션 필요

`YELLOW`/`RED`에는:
- 문제 요약
- 제안 레드라인 문구
- fallback
- 비즈니스 영향
- 에스컬레이션 경로

### 6. 레드라인 제안

각 이슈마다:
- 현재 문구
- 제안 문구
- 제안 사유(상대방 공유 가능)
- 협상 우선순위(`must-have`/`nice-to-have`)

### 7. 비즈니스 영향 요약

- 전체 리스크 프로필
- Top 3 이슈
- 협상 전략(선공/양보 항목)
- 일정상 고려사항

### 8. CLM 라우팅(연결 시)

리스크 수준에 맞는 승인 워크플로/라우팅 경로를 제안합니다.

## Output

- Contract Review Summary
- Key Findings
- Clause-by-Clause Analysis
- Negotiation Strategy
- Next Steps

## Notes

- 비영어 계약서는 언어 처리 전제 확인
- 50p+ 장문은 핵심 조항 우선 검토 후 전체 검토로 확장
- 법률 의사결정 전 전문가 검토 필요 문구를 항상 포함
