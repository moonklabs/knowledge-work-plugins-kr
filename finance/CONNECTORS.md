# 커넥터

## 도구 참조 방식

플러그인 파일에서는 `~~category`를 해당 카테고리에 사용자가 연결하는 도구의 플레이스홀더로 사용합니다. 예를 들어, `~~data warehouse`는 MCP 서버가 있는 Snowflake, BigQuery 또는 기타 웨어하우스를 의미할 수 있습니다.

플러그인은 **도구 비종속적(tool-agnostic)**입니다. 특정 제품이 아닌 카테고리(데이터 웨어하우스, 채팅, 프로젝트 트래커 등)를 기준으로 워크플로우를 기술합니다. `.mcp.json`에는 특정 MCP 서버가 사전 구성되어 있지만, 해당 카테고리의 모든 MCP 서버가 호환됩니다.

## 이 플러그인의 커넥터

| 카테고리 | 플레이스홀더 | 포함된 서버 | 기타 옵션 |
|----------|-------------|------------|-----------|
| 데이터 웨어하우스 | `~~data warehouse` | Snowflake\*, Databricks\*, BigQuery | Redshift, PostgreSQL |
| 이메일 | `~~email` | Microsoft 365 | — |
| 오피스 제품군 | `~~office suite` | Microsoft 365 | — |
| 채팅 | `~~chat` | Slack | Microsoft Teams |
| ERP / 회계 | `~~erp` | — (아직 지원되는 MCP 서버 없음) | NetSuite, SAP, QuickBooks, Xero |
| 분석 / BI | `~~analytics` | — (아직 지원되는 MCP 서버 없음) | Tableau, Looker, Power BI |

\* 플레이스홀더 — MCP URL 미구성
