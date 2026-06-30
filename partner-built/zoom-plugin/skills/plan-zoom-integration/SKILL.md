---
name: plan-zoom-integration
description: Zoom integration idea를 architecture, auth, delivery milestone이 포함된 implementation plan으로 바꿉니다. Practical build plan, phased delivery sequence, risk list, next-step recommendation이 필요할 때 사용합니다.
argument-hint: "<what you want to build>"
user-invocable: false
---

# /plan-zoom-integration

> If you see unfamiliar placeholders or need to check which tools are connected, see [CONNECTORS.md](../../CONNECTORS.md).

Create a practical build plan for a Zoom integration or app.

## Usage

```text
/plan-zoom-integration $ARGUMENTS
```

## Workflow

1. Capture the target user flow and success criteria.
2. Choose the correct Zoom surface and supporting services.
3. Define auth requirements, scopes, and account assumptions.
4. Break implementation into phases: prototype, core integration, reliability, and launch.
5. Call out hard risks early: OAuth setup, webhook verification, SDK environment limits, marketplace review, or MCP client constraints.
6. End with the smallest deliverable that proves the architecture.

## Output

- Architecture summary
- Zoom products and APIs required
- Auth and scope checklist
- Delivery phases
- Risks, open questions, and immediate next action

## Related Skills

- [start](../start/SKILL.md)
- [setup-zoom-oauth](../setup-zoom-oauth/SKILL.md)
- [build-zoom-meeting-app](../build-zoom-meeting-app/SKILL.md)
- [build-zoom-bot](../build-zoom-bot/SKILL.md)
