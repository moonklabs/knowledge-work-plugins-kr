# MCP 검색 및 연결

플러그인 커스터마이징 중 MCP를 찾고 연결하는 방법입니다.

## 사용 가능한 도구

### `search_mcp_registry`
MCP 디렉토리에서 사용 가능한 커넥터를 검색합니다.

**입력:** `{ "keywords": ["array", "of", "search", "terms"] }`

**출력:** 최대 10개의 결과가 반환되며, 각 결과에는 다음이 포함됩니다:
- `name`: MCP 표시 이름
- `description`: 한 줄 설명
- `tools`: MCP가 제공하는 도구 이름 목록
- `url`: MCP 엔드포인트 URL (`.mcp.json`에서 사용)
- `directoryUuid`: suggest_connectors에서 사용할 UUID
- `connected`: 사용자가 이 MCP에 연결되어 있는지 여부를 나타내는 불리언 값

### `suggest_connectors`
사용자가 MCP를 설치/연결할 수 있도록 연결 버튼을 표시합니다.

**입력:** `{ "directoryUuids": ["uuid1", "uuid2"] }`

**출력:** 각 MCP에 대한 연결 버튼이 있는 UI를 렌더링합니다

## 카테고리-키워드 매핑

| 카테고리 | 검색 키워드 |
|----------|-------------|
| `project-management` | `["asana", "jira", "linear", "monday", "tasks"]` |
| `software-coding` | `["github", "gitlab", "bitbucket", "code"]` |
| `chat` | `["slack", "teams", "discord"]` |
| `documents` | `["google docs", "notion", "confluence"]` |
| `calendar` | `["google calendar", "calendar"]` |
| `email` | `["gmail", "outlook", "email"]` |
| `design-graphics` | `["figma", "sketch", "design"]` |
| `analytics-bi` | `["datadog", "grafana", "analytics"]` |
| `crm` | `["salesforce", "hubspot", "crm"]` |
| `wiki-knowledge-base` | `["notion", "confluence", "outline", "wiki"]` |
| `data-warehouse` | `["bigquery", "snowflake", "redshift"]` |
| `conversation-intelligence` | `["gong", "chorus", "call recording"]` |

## 워크플로우

1. **커스터마이징 포인트 찾기**: `~~` 접두사가 붙은 값을 찾습니다(예: `~~Jira`)
2. **이전 단계 발견 사항 확인**: 어떤 도구를 사용하는지 이미 파악했는지 확인합니다
   - **예**: 해당 도구를 검색하여 `url`을 가져오고 5단계로 건너뜁니다
   - **아니오**: 3단계로 진행합니다
3. **검색**: 매핑된 키워드로 `search_mcp_registry`를 호출합니다
4. **선택지 제시 및 사용자 질문**: 모든 결과를 보여주고 어떤 것을 사용하는지 질문합니다
5. **필요 시 연결**: 연결되지 않은 경우 `suggest_connectors`를 호출합니다
6. **MCP 설정 업데이트**: 검색 결과의 `url`을 사용하여 설정을 추가합니다

## 플러그인 MCP 설정 업데이트

### 설정 파일 찾기

1. **`plugin.json`에서 `mcpServers` 필드를 확인합니다**:
   ```json
   {
     "name": "my-plugin",
     "mcpServers": "./config/servers.json"
   }
   ```
   해당 필드가 있으면 해당 경로의 파일을 편집합니다.

2. **`mcpServers` 필드가 없으면** 플러그인 루트의 `.mcp.json`을 사용합니다(기본값).

3. **`mcpServers`가 `.mcpb` 파일만 가리키는 경우**(번들 서버), 플러그인 루트에 새 `.mcp.json`을 생성합니다.

### 설정 파일 형식

래핑된 형식과 래핑되지 않은 형식 모두 지원됩니다:

```json
{
  "mcpServers": {
    "github": {
      "type": "http",
      "url": "https://api.githubcopilot.com/mcp/"
    }
  }
}
```

`search_mcp_registry` 결과의 `url` 필드를 사용합니다.

**참고:** 퍼스트파티 통합(Gmail, Google Calendar, Google Drive)은 사용자 수준에서 연결되므로 플러그인 `.mcp.json` 항목이 필요하지 않습니다.
