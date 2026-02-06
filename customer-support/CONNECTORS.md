# 커넥터

## 도구 레퍼런스 동작 방식

플러그인 파일은 `~~category`를 사용자가 연결한 해당 카테고리의 도구를 가리키는 플레이스홀더로 사용합니다. 예를 들어 `~~support platform`은 Intercom, Zendesk 또는 MCP 서버가 있는 다른 지원 도구를 의미할 수 있습니다.

플러그인은 **도구에 종속되지 않습니다**. 특정 제품명이 아니라 카테고리(지원 플랫폼, CRM, 채팅 등)로 워크플로를 설명합니다. `.mcp.json`에는 특정 MCP 서버가 미리 구성되어 있지만, 해당 카테고리의 어떤 MCP 서버라도 사용할 수 있습니다.

## 이 플러그인의 커넥터

| 카테고리 | 플레이스홀더 | 포함 서버 | 기타 옵션 |
|----------|-------------|-----------------|---------------|
| 채팅 | `~~chat` | Slack | Microsoft Teams |
| 이메일 | `~~email` | Microsoft 365 | — |
| 클라우드 스토리지 | `~~cloud storage` | Microsoft 365 | — |
| 지원 플랫폼 | `~~support platform` | Intercom | Zendesk, Freshdesk, HubSpot Service Hub |
| CRM | `~~CRM` | HubSpot | Salesforce, Pipedrive |
| 지식베이스 | `~~knowledge base` | Guru, Notion | Confluence, Help Scout |
| 프로젝트 트래커 | `~~project tracker` | Atlassian (Jira/Confluence) | Linear, Asana |
