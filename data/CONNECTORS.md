# 커넥터

## 도구 참조 방식

플러그인 파일은 해당 카테고리에서 사용자가 연결한 도구의 플레이스홀더로 `~~category`를 사용합니다. 예를 들어, `~~data warehouse`는 MCP 서버가 있는 Snowflake, BigQuery 또는 기타 웨어하우스를 의미할 수 있습니다.

플러그인은 **도구에 구애받지 않습니다** -- 특정 제품이 아닌 카테고리(데이터 웨어하우스, 노트북, 제품 분석 등)를 기준으로 워크플로우를 설명합니다. `.mcp.json`은 특정 MCP 서버를 사전 구성하지만, 해당 카테고리의 모든 MCP 서버가 작동합니다.

## 이 플러그인의 커넥터

| 카테고리 | 플레이스홀더 | 포함된 서버 | 기타 옵션 |
|----------|-------------|------------|----------|
| 데이터 웨어하우스 | `~~data warehouse` | Snowflake\*, Databricks\*, BigQuery | Redshift, PostgreSQL, MySQL |
| 노트북 | `~~notebook` | Hex | Jupyter, Deepnote, Observable |
| 제품 분석 | `~~product analytics` | Amplitude | Mixpanel, Heap |
| 프로젝트 트래커 | `~~project tracker` | Atlassian (Jira/Confluence) | Linear, Asana |

\* 플레이스홀더 -- MCP URL이 아직 구성되지 않음
