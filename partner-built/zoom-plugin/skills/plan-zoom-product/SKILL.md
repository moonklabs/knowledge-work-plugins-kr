---
name: plan-zoom-product
description: >-
  사용 사례에 맞는 Zoom 빌드 표면을 선택하고 트레이드오프를 설명합니다. 구체적인 제품 아이디어나 통합 목표에 대해 REST API, Webhooks,
  WebSockets, Meeting SDK, Video SDK, Zoom Apps SDK, Phone, Contact Center, MCP 중에서 결정할 때 사용합니다.
argument-hint: "<product idea, app type, or integration goal>"
user-invocable: false
---

# /plan-zoom-product

> If you see unfamiliar placeholders or need to check which tools are connected, see [CONNECTORS.md](../../CONNECTORS.md).

Choose between Zoom REST API, Webhooks, WebSockets, Meeting SDK, Video SDK, Zoom Apps SDK, Phone, Contact Center, or MCP for a specific use case.

## Usage

```text
/plan-zoom-product $ARGUMENTS
```

## Workflow

1. Identify the user's actual goal.
2. Classify whether the problem is automation, embedded meetings, custom video, in-client app behavior, event delivery, AI tooling, or support/phone/contact-center work.
3. If the request is ambiguous, ask one short clarifier before locking the recommendation.
4. Recommend the primary Zoom surface and list the minimum supporting pieces.
5. Explain why the rejected alternatives are worse for this case.
6. End with a concrete next-step plan.

## Output

- Recommended Zoom surface
- Supporting components required
- Key tradeoffs and constraints
- Suggested implementation sequence
- Relevant skill links for the next step

## Related Skills

- [start](../start/SKILL.md)
- [choose-zoom-approach](../choose-zoom-approach/SKILL.md)
- [design-mcp-workflow](../design-mcp-workflow/SKILL.md)
