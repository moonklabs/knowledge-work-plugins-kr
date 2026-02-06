# 커넥터

## 도구 레퍼런스 동작 방식

플러그인 파일은 `~~category`를 사용자가 연결한 해당 카테고리의 도구를 가리키는 플레이스홀더로 사용합니다. 예를 들어 `~~data warehouse`는 Snowflake, BigQuery 또는 MCP 서버가 있는 다른 웨어하우스를 의미할 수 있습니다.

플러그인은 **도구에 종속되지 않습니다**. 특정 제품명이 아니라 카테고리(데이터 웨어하우스, 노트북, 제품 분석 등)로 워크플로를 설명합니다. `.mcp.json`에는 특정 MCP 서버가 미리 구성되어 있지만, 해당 카테고리의 어떤 MCP 서버라도 사용할 수 있습니다.

## 이 플러그인의 커넥터

| 카테고리 | 플레이스홀더 | 포함 서버 | 기타 옵션 |
|----------|-------------|-----------------|---------------|
| 데이터 웨어하우스 | `~~data warehouse` | Snowflake\*, Databricks\*, BigQuery | Redshift, PostgreSQL, MySQL |
| 노트북 | `~~notebook` | Hex | Jupyter, Deepnote, Observable |
| 제품 분석 | `~~product analytics` | Amplitude | Mixpanel, Heap |
| 프로젝트 트래커 | `~~project tracker` | Atlassian (Jira/Confluence) | Linear, Asana |

\* 플레이스홀더 — MCP URL이 아직 설정되지 않았습니다
