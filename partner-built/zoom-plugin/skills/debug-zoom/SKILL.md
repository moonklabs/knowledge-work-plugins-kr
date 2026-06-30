---
name: debug-zoom
description: >-
  실패 지점을 격리하고 적절한 Zoom 참조 스킬로 라우팅해 깨진 Zoom 통합을 디버그합니다. Auth, API, 웹훅, SDK, MCP 동작이 실패하고 우선순위가 있는
  가설 목록과 검증 단계가 필요할 때 사용합니다.
argument-hint: "<symptoms, error, or failing flow>"
---

# /debug-zoom

> If you see unfamiliar placeholders or need to check which tools are connected, see [CONNECTORS.md](../../CONNECTORS.md).

Debug Zoom auth, API, webhook, SDK, or MCP issues without wandering through the entire docs set.

## Usage

```text
/debug-zoom $ARGUMENTS
```

## Workflow

1. Identify the failing layer: auth, API request, webhook, SDK init, media/session behavior, or MCP transport.
2. Ask for the minimum missing evidence: exact error, platform, request/response, event payload, or code path.
3. Produce 2-4 plausible causes ranked by likelihood.
4. Route to the most relevant deep references in `skills/`.
5. Give a short verification plan so the user can confirm the fix.

## Output

- Most likely failure layer
- Ranked hypotheses
- Targeted fix steps
- Verification checklist
- Relevant skill links

## Related Skills

- [debug-zoom-integration](../debug-zoom-integration/SKILL.md)
- [setup-zoom-oauth](../setup-zoom-oauth/SKILL.md)
- [design-mcp-workflow](../design-mcp-workflow/SKILL.md)
