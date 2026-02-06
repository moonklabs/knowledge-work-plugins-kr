---
description: 수신 NDA를 신속 분류(표준 승인/법무 검토/전면 검토)
argument-hint: "<NDA file or text>"
---

# /triage-nda -- NDA Pre-Screening

> 낯선 플레이스홀더가 보이거나 어떤 도구가 연결되어 있는지 확인하려면 [CONNECTORS.md](../CONNECTORS.md)를 보세요.

수신 NDA를 기준표로 빠르게 스크리닝하고 라우팅 등급을 부여합니다.

## Invocation

```bash
/triage-nda
```

## Workflow

### 1. NDA 입력

- 파일 업로드
- URL
- 텍스트 붙여넣기

### 2. NDA 플레이북 로드

`legal.local.md`에서 기준을 찾습니다.
기준이 없으면 시장 표준 기본값 적용:
- 상호 의무 선호
- 기간 2-3년(영업비밀 최대 5년)
- 표준 carveout 필수
- 비경쟁/비권유 금지
- 과도한 residuals 제한
- 합리적 관할

### 3. 퀵 스크린

점검 항목:
- 상호/편무 의무
- 비밀정보 정의 범위
- 기간
- carveout 완비 여부
- 허용된 공유 대상
- 반환/파기 의무
- residuals
- 비권유/비경쟁 포함 여부
- 금지명령 조항 균형성
- 준거법
- 양도
- NDA 외 이례 조항

### 4. 분류

- `GREEN`: 표준 승인 가능
- `YELLOW`: 법무 리뷰 필요
- `RED`: 중대 이슈, 전면 검토 필요

### 5. 트리아지 리포트 생성

포함:
- 분류 결과
- 당사자/유형/기간/준거법
- 항목별 PASS/FLAG/FAIL
- 이슈별 리스크와 수정 제안
- 다음 조치

### 6. 라우팅 제안

- GREEN: 표준 권한으로 진행
- YELLOW: 지정 이슈 중심 검토 라우팅
- RED: 전면 검토 및 표준 NDA 반제안 권고

## Notes

- NDA가 아닌 상업 조건이 섞인 문서면 즉시 RED
- MSA 일부로 포함된 NDA 조항은 본계약 맥락 반영 필요
- 본 커맨드는 스크리닝 도구이며 확정 판단은 법무 검토 필요
