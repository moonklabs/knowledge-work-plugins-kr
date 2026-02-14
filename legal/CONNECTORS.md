# 커넥터

## 도구 참조 방식

플러그인 파일은 해당 카테고리에서 사용자가 연결하는 도구에 대한 플레이스홀더로 `~~category`를 사용합니다. 예를 들어, `~~cloud storage`는 MCP 서버가 있는 Box, Egnyte 또는 기타 스토리지 제공자를 의미할 수 있습니다.

플러그인은 **도구에 종속되지 않으며**, 특정 제품이 아닌 카테고리(클라우드 스토리지, 채팅, 오피스 제품군 등)로 워크플로우를 기술합니다. `.mcp.json`에 특정 MCP 서버가 사전 구성되어 있지만, 해당 카테고리의 모든 MCP 서버가 호환됩니다.

## 이 플러그인의 커넥터

| 카테고리 | 플레이스홀더 | 포함된 서버 | 기타 옵션 |
|----------|-------------|------------|----------|
| 채팅 | `~~chat` | Slack | Microsoft Teams |
| 클라우드 스토리지 | `~~cloud storage` | Box, Egnyte | Dropbox, SharePoint, Google Drive |
| CLM | `~~CLM` | — | Ironclad, Agiloft |
| CRM | `~~CRM` | — | Salesforce, HubSpot |
| 전자서명 | `~~e-signature` | — | DocuSign, Adobe Sign |
| 오피스 제품군 | `~~office suite` | Microsoft 365 | Google Workspace |
| 프로젝트 트래커 | `~~project tracker` | Atlassian (Jira/Confluence) | Linear, Asana |
