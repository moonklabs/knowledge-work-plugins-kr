# 커넥터

## 도구 참조 방식

플러그인 파일은 `~~category`를 해당 카테고리에서 사용자가 연결한 도구의 플레이스홀더로 사용합니다. 예를 들어 `~~CRM`은 Salesforce, HubSpot 또는 MCP 서버가 있는 다른 CRM을 의미할 수 있습니다.

플러그인은 **도구에 구애받지 않습니다** — 특정 제품이 아닌 카테고리(CRM, 채팅, 이메일 등)로 워크플로우를 설명합니다. `.mcp.json`에 특정 MCP 서버가 사전 구성되어 있지만, 해당 카테고리의 어떤 MCP 서버든 사용할 수 있습니다.

## 이 플러그인의 커넥터

| 카테고리 | 플레이스홀더 | 포함 서버 | 기타 옵션 |
|----------|-------------|-----------|-----------|
| 캘린더 | `~~calendar` | Microsoft 365 | Google Calendar |
| 채팅 | `~~chat` | Slack | Microsoft Teams |
| CRM | `~~CRM` | HubSpot, Close | Salesforce, Pipedrive, Copper |
| 데이터 보강 | `~~data enrichment` | Clay, ZoomInfo | Apollo, Clearbit, Lusha |
| 이메일 | `~~email` | Microsoft 365 | Gmail |
| 지식 베이스 | `~~knowledge base` | Notion | Confluence, Guru |
| 미팅 녹취 | `~~conversation intelligence` | Fireflies | Gong, Chorus, Otter.ai |
| 프로젝트 트래커 | `~~project tracker` | Atlassian (Jira/Confluence) | Linear, Asana |
