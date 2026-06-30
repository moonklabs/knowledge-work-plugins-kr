# Zoom 플러그인

Zoom 통합을 계획하고 구축하며 디버그하기 위한 Claude 플러그인입니다. 사용자가 전체 문서 트리를 먼저 읽지 않아도 적절한 Zoom 구현 영역을 고르고, 구현을 구체화하고, 실패를 디버그하며, 필요한 Zoom 참조로 라우팅할 수 있게 돕습니다.

## 설치

이 디렉터리를 로컬 Claude 플러그인으로 설치합니다. 플러그인 매니페스트는 [`.claude-plugin/plugin.json`](.claude-plugin/plugin.json)에 있고, 함께 제공되는 Zoom MCP 커넥터는 [`.mcp.json`](.mcp.json)에 정의되어 있습니다.

함께 제공되는 MCP 서버를 사용하기 전에 Claude가 사용할 Zoom 구현 영역용 bearer token을 내보내세요.

```bash
export ZOOM_MCP_ACCESS_TOKEN="your_zoom_user_oauth_access_token"
export ZOOM_DOCS_MCP_ACCESS_TOKEN="your_zoom_docs_mcp_access_token"
export ZOOM_WHITEBOARD_MCP_ACCESS_TOKEN="your_zoom_user_oauth_access_token"
```

## 슬래시 워크플로

`skills/` 아래 스킬로 구현된 명시적 슬래시 워크플로입니다.

| 워크플로 | 설명 |
|---|---|
| [`/start`](skills/start/SKILL.md) | Zoom 앱 아이디어에서 시작해 적절한 제품과 빌드 경로로 라우팅합니다 |
| [`/setup-zoom-oauth`](skills/setup-zoom-oauth/SKILL.md) | Zoom app용 인증 모델, scope, redirect flow를 선택합니다 |
| [`/build-zoom-meeting-app`](skills/build-zoom-meeting-app/SKILL.md) | 임베드형 또는 관리형 Zoom 미팅 흐름을 구축합니다 |
| [`/build-zoom-bot`](skills/build-zoom-bot/SKILL.md) | 봇, 레코더, 실시간 미팅 처리기를 구축합니다 |
| [`/debug-zoom`](skills/debug-zoom/SKILL.md) | 깨진 Zoom 통합을 분류하고 실패 계층을 격리합니다 |
| [`/setup-zoom-mcp`](skills/setup-zoom-mcp/SKILL.md) | Zoom MCP가 적합한 시점을 판단하고 안전한 Claude 워크플로를 설정합니다 |
| [`/build-zoom-rest-api-app`](skills/rest-api/SKILL.md) | Zoom REST endpoint, scope, resource pattern으로 라우팅합니다 |
| [`/build-zoom-meeting-sdk-app`](skills/meeting-sdk/SKILL.md) | 임베드형 Zoom 미팅 구현 세부사항으로 라우팅합니다 |
| [`/build-zoom-video-sdk-app`](skills/video-sdk/SKILL.md) | 커스텀 video-session 구현 세부사항으로 라우팅합니다 |
| [`/setup-zoom-webhooks`](skills/webhooks/SKILL.md) | Zoom webhook subscription, signature verification, handler를 설정합니다 |
| [`/setup-zoom-websockets`](skills/websockets/SKILL.md) | Webhook보다 적합할 때 Zoom WebSocket event delivery를 설정합니다 |
| [`/build-zoom-team-chat-app`](skills/team-chat/SKILL.md) | Team Chat 사용자 또는 챗봇 통합을 구축합니다 |
| [`/build-zoom-phone-integration`](skills/phone/SKILL.md) | Smart Embed, API, event 기반 Zoom Phone 통합을 구축합니다 |
| [`/build-zoom-contact-center-app`](skills/contact-center/SKILL.md) | Contact Center 앱, 웹, 네이티브 통합을 구축합니다 |
| [`/build-zoom-virtual-agent`](skills/virtual-agent/SKILL.md) | Virtual Agent 웹 또는 모바일 래퍼 통합을 구축합니다 |

## 내부 라우팅 스킬

이 항목들은 자동 라우팅 도우미로 플러그인에 남아 있지만 공개 슬래시 명령 표면에는 포함되지 않습니다.

- [`start`](skills/start/SKILL.md)
- [`plan-zoom-product`](skills/plan-zoom-product/SKILL.md)
- [`plan-zoom-integration`](skills/plan-zoom-integration/SKILL.md)
- [`choose-zoom-approach`](skills/choose-zoom-approach/SKILL.md)
- [`design-mcp-workflow`](skills/design-mcp-workflow/SKILL.md)
- [`debug-zoom-integration`](skills/debug-zoom-integration/SKILL.md)

## 심층 참조

플러그인은 원래의 Zoom 제품별 참조 라이브러리도 `skills/` 아래에 유지합니다. 이는 기본 진입 표면이 아니라 보조 참조입니다.

- [`skills/general/`](skills/general/)
- [`skills/rest-api/`](skills/rest-api/)
- [`skills/meeting-sdk/`](skills/meeting-sdk/)
- [`skills/video-sdk/`](skills/video-sdk/)
- [`skills/webhooks/`](skills/webhooks/)
- [`skills/websockets/`](skills/websockets/)
- [`skills/oauth/`](skills/oauth/)
- [`skills/zoom-mcp/`](skills/zoom-mcp/)

## 예시 워크플로

### Zoom 앱 아이디어에서 시작

```text
/start 통화에 참여해 액션 아이템을 추출하고 요약을 저장하는 내부 미팅 비서를 만들어줘
```

### 새 앱 계획

```text
/start 고객이 우리 제품에서 Zoom 미팅을 예약하고 참여할 수 있는 React app을 만들어줘
```

### 깨진 웹훅 디버깅

```text
/debug-zoom Zoom webhook 서명 검증이 로컬에서는 되는데 production에서 실패해
```

### MCP 흐름 설계

```text
/setup-zoom-mcp Claude가 미팅을 검색하고 녹화 리소스를 가져온 뒤 후속 문서를 만들게 하고 싶어
```

## 커넥터

[CONNECTORS.md](CONNECTORS.md)를 참고하세요. 플러그인은 함께 제공되는 스킬만으로 단독 동작하며, Claude가 [`.mcp.json`](.mcp.json)의 bundled Zoom MCP 서버를 사용할 수 있으면 더 강력해집니다.

## 크로스 플랫폼 참고

이 저장소는 우선 Claude 플러그인으로 패키징되어 있지만, 저장소 수준 discovery file을 사용하는 agent ecosystem을 위해 [AGENTS.md](AGENTS.md)도 포함합니다. 재사용 가능한 핵심은 `skills/` 트리와 그 안의 `SKILL.md` 파일입니다.

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
│   └── ... 기존 Zoom 참조 스킬
└── README.md
```

## 라이선스

MIT
