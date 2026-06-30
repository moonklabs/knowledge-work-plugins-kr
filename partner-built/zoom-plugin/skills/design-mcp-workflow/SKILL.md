---
name: design-mcp-workflow
description: Claude를 위한 Zoom MCP 워크플로를 설계합니다. Zoom MCP가 작업에 맞는지 판단하거나, tool-기반 AI 워크플로를 계획하거나, MCP responsibility와 REST API responsibility를 분리할 때 사용합니다.
user-invocable: false
---

# Design MCP Workflow

Use this skill when the user wants Claude or another MCP-capable client to interact with Zoom via tool calls instead of only deterministic API code.

## Covers

- MCP fit assessment
- REST API vs MCP boundaries
- Hybrid architectures
- Connector expectations
- Whiteboard-specific MCP routing

## Workflow

1. Decide whether the problem is agentic tooling, deterministic automation, or both.
2. Route MCP-only tasks to [zoom-mcp](../zoom-mcp/SKILL.md).
3. Route hybrid tasks to both [zoom-mcp](../zoom-mcp/SKILL.md) and [rest-api](../rest-api/SKILL.md).
4. If Whiteboard is central, route to [zoom-mcp/whiteboard](../zoom-mcp/whiteboard/SKILL.md).
5. Call out transport, auth, and client capability assumptions explicitly.

## Common Mistakes

- Using MCP for deterministic backend jobs that should stay in REST
- Treating MCP as a replacement for all API design
- Ignoring client transport support and auth requirements
