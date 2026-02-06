---
description: 콘텐츠를 브랜드 보이스, 스타일 가이드, 메시지 기둥 기준으로 검토
argument-hint: "<content to review>"
---

# Brand Review

> 낯선 플레이스홀더가 보이거나 어떤 도구가 연결되어 있는지 확인하려면 [CONNECTORS.md](../CONNECTORS.md)를 보세요.

마케팅 콘텐츠를 브랜드 기준으로 검토하고 이탈 항목과 수정안을 제공합니다.

## Trigger

```bash
/brand-review <content>
```

## Workflow

1. 브랜드 기준 수집: 로컬 설정(보이스, 금지 표현, 핵심 메시지, 용어집)을 로드합니다.
2. 콘텐츠 분석: 톤, 구조, 메시지 일관성, CTA, 채널 적합성을 점검합니다.
3. 이슈 분류: `Must fix`, `Should fix`, `Optional`로 분류합니다.
4. 리라이트 제안: 문장 단위 대체 문구와 이유를 제공합니다.
5. 최종판 제안: 수정 반영 버전과 짧은 변경 요약을 제공합니다.

## Output

- Overall Fit Score
- Findings by category(voice, clarity, messaging, compliance)
- Priority fixes
- Revised draft
