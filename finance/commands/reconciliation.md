---
description: GL 잔액을 서브레저/은행/외부 잔액과 대사
argument-hint: "<account> [period]"
---

# Account Reconciliation

> 낯선 플레이스홀더가 보이거나 어떤 도구가 연결되어 있는지 확인하려면 [CONNECTORS.md](../CONNECTORS.md)를 보세요.

**중요**: 이 커맨드는 대사 워크플로를 보조하지만 재무 자문을 제공하지 않습니다. 사인오프 전 전문가 검토가 필요합니다.

## Usage

```bash
/recon <account> <period>
```

### Arguments

- `account`: `cash`, `ar`, `ap`, `fixed-assets`, `intercompany`, `prepaid`, `accrued-liabilities`, 또는 계정코드
- `period`: `2024-12`, `2024-12-31` 등

## Workflow

### 1. 양측 잔액 수집

연결 소스가 있으면:
- 기간말 GL 잔액
- 비교측 잔액(서브레저/은행/외부확인)
- 전기 미결항목

연결이 없으면 사용자에게 요청:
1. GL 잔액
2. 비교측 잔액
3. 전기 미결항목(선택)

### 2. 차이 계산

`GL Balance - Other Balance`를 계산합니다.

### 3. 조정항목 식별

- 타이밍 차이: 미결제수표, 미입금, 반제 대기 등
- 영구 차이: 오분개, 누락분개, 수수료 미반영 등
- 전기 이월: 해결/미해결 구분 및 에이징

### 4. 대사 워크페이퍼 생성

다음을 포함:
- 계정/기간/작성자
- GL 기준 잔액
- 가산/차감 조정항목
- 조정 후 잔액
- 비교측 잔액 및 조정
- 최종 차이(`$0.00` 목표)

### 5. 조정항목 상세

`Description`, `Amount`, `Category`, `Age`, `Status`, `Action Required` 표로 제시합니다.

### 6. 에스컬레이션 플래그

- 장기 미해결(30/60/90일)
- 중요 금액 초과
- 전기 대비 증가 추세
- 원인 미확인 차이

### 7. 출력

1. 대사 워크페이퍼
2. 조정항목 목록/에이징
3. 필요 수정분개
4. 후속 액션
5. 전기 대비
6. 작성자/검토자 사인오프 섹션
