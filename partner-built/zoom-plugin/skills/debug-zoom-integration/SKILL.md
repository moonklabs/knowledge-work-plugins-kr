---
name: debug-zoom-integration
description: >-
  깨진 Zoom 구현을 빠르게 디버그합니다. Auth, 웹훅, SDK 조인, MCP 전송, 실시간 미디어 워크플로가 실패했고 수정안을 내기 전에 계층을 먼저 격리해야 할 때
  사용합니다.
user-invocable: false
---

# Debug Zoom Integration

Use this skill when the user already built something and it is failing.

## Triage Order

1. Auth and app configuration
2. Request construction or event verification
3. SDK initialization or platform mismatch
4. Media/session behavior
5. MCP transport and capability assumptions

## Evidence To Request

- Exact error text
- Platform and SDK/runtime
- Relevant request or payload sample
- What worked versus what failed
- Whether the issue is reproducible or intermittent

## Reference Routing

- [oauth](../oauth/SKILL.md)
- [rest-api](../rest-api/SKILL.md)
- [webhooks](../webhooks/SKILL.md)
- [meeting-sdk](../meeting-sdk/SKILL.md)
- [video-sdk](../video-sdk/SKILL.md)
- [rtms](../rtms/SKILL.md)
- [zoom-mcp](../zoom-mcp/SKILL.md)

## Output

- Most likely failing layer
- Ranked hypotheses
- Short fix plan
- Verification steps
