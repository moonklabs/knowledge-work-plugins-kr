---
description: 지원 티켓/고객 이슈 분류 및 우선순위 지정
argument-hint: "<ticket or issue description>"
---

# Triage

> 낯선 플레이스홀더가 보이거나 어떤 도구가 연결되어 있는지 확인하려면 [CONNECTORS.md](../CONNECTORS.md)를 보세요.

수신 이슈를 카테고리/우선순위/라우팅으로 구조화합니다.

## Usage

```bash
/triage <issue>
```

## Workflow

1. 이슈 요약 및 영향 범위 파악
2. 카테고리/심각도(P1-P4) 분류
3. 라우팅 대상(Tier, Engineering, Product) 추천
4. 초기 응답 초안 생성
5. 추적 필드(재현성, 계정 영향, SLA) 정리

## Output

- Triage classification
- Priority rationale
- Routing recommendation
- Initial response draft
