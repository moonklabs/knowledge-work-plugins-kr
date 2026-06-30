# Zoom 플러그인

Zoom integration을 계획하고 구축하며 debugging하기 위한 Claude plugin입니다. 사용자가 전체 doc tree를 먼저 읽지 않아도 적절한 Zoom surface를 고르고, implementation을 구체화하고, failure를 debug하며, 필요한 Zoom reference로 라우팅할 수 있게 돕습니다.

## 설치

이 directory를 local Claude plugin으로 설치합니다. Plugin manifest는 [`.claude-plugin/plugin.json`](.claude-plugin/plugin.json)에 있고, bundled Zoom MCP connector는 [`.mcp.json`](.mcp.json)에 정의되어 있습니다.

Bundled MCP server를 사용하기 전에 Claude가 사용할 Zoom surface용 bearer token을 export하세요.

```bash
export ZOOM_MCP_ACCESS_TOKEN="your_zoom_user_oauth_access_token"
export ZOOM_DOCS_MCP_ACCESS_TOKEN="your_zoom_docs_mcp_access_token"
export ZOOM_WHITEBOARD_MCP_ACCESS_TOKEN="your_zoom_user_oauth_access_token"
```

## Slash workflow

`skills/` 아래 skill로 구현된 explicit slash workflow입니다.

| Workflow | 설명 |
|---|---|
| [`/start`](skills/start/SKILL.md) | Zoom app idea에서 시작해 적절한 product와 build path로 routing합니다 |
| [`/setup-zoom-oauth`](skills/setup-zoom-oauth/SKILL.md) | Zoom app용 auth model, scope, redirect flow를 선택합니다 |
| [`/build-zoom-meeting-app`](skills/build-zoom-meeting-app/SKILL.md) | Embedded 또는 managed Zoom meeting flow를 구축합니다 |
| [`/build-zoom-bot`](skills/build-zoom-bot/SKILL.md) | Bot, recorder, real-time meeting processor를 구축합니다 |
| [`/debug-zoom`](skills/debug-zoom/SKILL.md) | Broken Zoom integration을 triage하고 failing layer를 격리합니다 |
| [`/setup-zoom-mcp`](skills/setup-zoom-mcp/SKILL.md) | Zoom MCP가 적합한 시점을 판단하고 안전한 Claude workflow를 설정합니다 |
| [`/build-zoom-rest-api-app`](skills/rest-api/SKILL.md) | Zoom REST endpoint, scope, resource pattern으로 routing합니다 |
| [`/build-zoom-meeting-sdk-app`](skills/meeting-sdk/SKILL.md) | Embedded Zoom meeting implementation detail로 routing합니다 |
| [`/build-zoom-video-sdk-app`](skills/video-sdk/SKILL.md) | Custom video-session implementation detail로 routing합니다 |
| [`/setup-zoom-webhooks`](skills/webhooks/SKILL.md) | Zoom webhook subscription, signature verification, handler를 설정합니다 |
| [`/setup-zoom-websockets`](skills/websockets/SKILL.md) | Webhook보다 적합할 때 Zoom WebSocket event delivery를 설정합니다 |
| [`/build-zoom-team-chat-app`](skills/team-chat/SKILL.md) | Team Chat user 또는 chatbot integration을 구축합니다 |
| [`/build-zoom-phone-integration`](skills/phone/SKILL.md) | Smart Embed, API, event 기반 Zoom Phone integration을 구축합니다 |
| [`/build-zoom-contact-center-app`](skills/contact-center/SKILL.md) | Contact Center app, web, native integration을 구축합니다 |
| [`/build-zoom-virtual-agent`](skills/virtual-agent/SKILL.md) | Virtual Agent web 또는 mobile wrapper integration을 구축합니다 |

## 내부 routing skill

이 항목들은 automatic routing helper로 plugin에 남아 있지만 public slash-command surface에는 포함되지 않습니다.

- [`start`](skills/start/SKILL.md)
- [`plan-zoom-product`](skills/plan-zoom-product/SKILL.md)
- [`plan-zoom-integration`](skills/plan-zoom-integration/SKILL.md)
- [`choose-zoom-approach`](skills/choose-zoom-approach/SKILL.md)
- [`design-mcp-workflow`](skills/design-mcp-workflow/SKILL.md)
- [`debug-zoom-integration`](skills/debug-zoom-integration/SKILL.md)

## Deep reference

Plugin은 original Zoom product-specific reference library도 `skills/` 아래에 유지합니다. 이는 primary entry surface가 아니라 supporting reference입니다.

- [`skills/general/`](skills/general/)
- [`skills/rest-api/`](skills/rest-api/)
- [`skills/meeting-sdk/`](skills/meeting-sdk/)
- [`skills/video-sdk/`](skills/video-sdk/)
- [`skills/webhooks/`](skills/webhooks/)
- [`skills/websockets/`](skills/websockets/)
- [`skills/oauth/`](skills/oauth/)
- [`skills/zoom-mcp/`](skills/zoom-mcp/)

## 예시 워크플로

### Zoom app idea에서 시작

```text
/start Build an internal meeting assistant that joins calls, extracts action items, and stores summaries
```

### 새 app 계획

```text
/start Build a React app that lets customers schedule and join Zoom meetings from our product
```

### Broken webhook debugging

```text
/debug-zoom My Zoom webhook signature verification fails in production but not locally
```

### MCP flow 설계

```text
/setup-zoom-mcp I want Claude to search meetings, pull recording resources, and create follow-up docs
```

## Connector

[CONNECTORS.md](CONNECTORS.md)를 참고하세요. Plugin은 bundled skill만으로 standalone 동작하며, Claude가 [`.mcp.json`](.mcp.json)의 bundled Zoom MCP server를 사용할 수 있으면 더 강력해집니다.

## Cross-platform notes

이 repo는 우선 Claude plugin으로 package되어 있지만, repo-level discovery file을 사용하는 agent ecosystem을 위해 [AGENTS.md](AGENTS.md)도 포함합니다. Reusable core는 `skills/` tree와 그 안의 `SKILL.md` files입니다.

## 구조

```text
Zoom Plugin/
├── .claude-plugin/plugin.json
├── .mcp.json
├── CONNECTORS.md
├── skills/
│   ├── plan-zoom-product/
│   ├── plan-zoom-integration/
│   ├── debug-zoom/
│   ├── setup-zoom-mcp/
│   ├── start/
│   ├── choose-zoom-approach/
│   ├── setup-zoom-oauth/
│   ├── build-zoom-meeting-app/
│   ├── build-zoom-bot/
│   ├── design-mcp-workflow/
│   ├── debug-zoom-integration/
│   └── ... existing Zoom reference skills
└── README.md
```

## 라이선스

MIT
