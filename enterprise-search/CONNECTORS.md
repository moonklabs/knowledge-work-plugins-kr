# 커넥터

## 도구 레퍼런스 동작 방식

플러그인 파일은 `~~category`를 사용자가 연결한 해당 카테고리의 도구를 가리키는 플레이스홀더로 사용합니다. 예를 들어 `~~chat`은 Slack, Microsoft Teams 또는 MCP 서버가 있는 다른 채팅 도구를 의미할 수 있습니다.

플러그인은 **도구에 종속되지 않습니다**. 특정 제품명이 아니라 카테고리(채팅, 이메일, 클라우드 스토리지 등)로 워크플로를 설명합니다. `.mcp.json`에는 특정 MCP 서버가 미리 구성되어 있지만, 해당 카테고리의 어떤 MCP 서버라도 사용할 수 있습니다.

이 플러그인은 검색 출력에서 소스 라벨로 `~~category` 참조를 광범위하게 사용합니다(예: `~~chat:`, `~~email:`). 이는 의도된 동작이며, 연결된 도구에 따라 동적으로 해석되는 카테고리 마커입니다.

## 이 플러그인의 커넥터

| 카테고리 | 플레이스홀더 | 포함 서버 | 기타 옵션 |
|----------|-------------|-----------------|---------------|
| 채팅 | `~~chat` | Slack | Microsoft Teams, Discord |
| 이메일 | `~~email` | Microsoft 365 | — |
| 클라우드 스토리지 | `~~cloud storage` | Microsoft 365 | Dropbox |
| 지식베이스 | `~~knowledge base` | Notion, Guru | Confluence, Slite |
| 프로젝트 트래커 | `~~project tracker` | Atlassian (Jira/Confluence), Asana | Linear, monday.com |
| CRM | `~~CRM` | *(사전 구성 없음)* | Salesforce, HubSpot |
| 오피스 제품군 | `~~office suite` | Microsoft 365 | Google Workspace |
