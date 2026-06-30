# Claude Code 및 Cowork용 Apollo 플러그인

[Apollo.io](https://www.apollo.io/)로 prospect를 찾고, lead를 enrich하고, outreach sequence를 load합니다. Apollo MCP Server 기반이며 **one-click integration**을 지원합니다.

---

## One-click MCP server 통합

이 plugin은 설치 시 **Apollo MCP Server를 자동으로 구성**합니다. Manual server setup이나 config file 수정 없이 plugin을 설치하고 Apollo Account로 인증하면 됩니다.

---

## 주요 스킬

이 plugin은 여러 Apollo API를 complete workflow로 연결하는 high-value skill을 제공합니다.

| 스킬 | 수행 내용 |
|---|---|
| `/apollo:enrich-lead` | Name, LinkedIn URL, email을 넣으면 email, phone, company intel, next action이 포함된 full contact card를 받습니다 |
| `/apollo:prospect` | ICP를 자연어로 설명하면 enriched decision-maker lead의 ranked table을 받습니다 |
| `/apollo:sequence-load` | Lead를 찾고 enrich한 뒤 outreach sequence에 bulk-load합니다. Dedup과 enrollment를 처리합니다 |

### `/apollo:enrich-lead`

Name, company, LinkedIn URL, email을 넣으면 email, phone, title, location, company detail, suggested next action이 포함된 complete contact card를 받습니다. "CEO of Figma" 같은 fuzzy lookup을 처리하고 exact match가 실패하면 search로 fallback합니다.

### `/apollo:prospect`

ICP를 자연어로 설명하세요. Pipeline이 matching company를 검색하고 firmographic data를 bulk-enrich하며 decision maker를 찾고 bulk enrichment로 contact info를 reveal한 뒤 ICP fit score가 포함된 ranked lead table을 반환합니다.

### `/apollo:sequence-load`

Targeting criteria에 맞는 contact를 찾고 enrich한 뒤 deduplication과 함께 contact로 생성하고 existing Apollo sequence에 bulk-add합니다. Enrollment 전에 candidate를 preview하고 이후 full summary를 보여줍니다.

---

## 설치

### Cowork

아래 link를 클릭해 한 번에 설치합니다.

[Install in Cowork](https://claude.ai/desktop/customize/plugins/new?marketplace=apolloio/apollo-mcp-plugin&plugin=apollo)

MCP server가 올바르게 시작되도록 Cowork를 다시 시작하세요.

### Claude Code

#### 1. 이 plugin marketplace 추가

Claude Code에서 실행합니다.

```
/plugin marketplace add apolloio/apollo-mcp-plugin
```

#### 2. Plugin 설치

```
/plugin install apollo@apollo-plugin-marketplace
```

#### 3. Claude Code 재시작

이렇게 하면 MCP server가 올바르게 시작됩니다.

---

## 인증

The Apollo MCP Server supports **OAuth**:

1. 설치 후 Claude Code 또는 Cowork에서 `/mcp`를 실행합니다
2. **Apollo** server를 선택하고 **Authenticate**를 클릭합니다
3. Browser에서 Apollo.io login을 완료합니다
4. 완료되면 모든 command를 사용할 준비가 됩니다

---

## Apollo credit

Some operations consume [Apollo credits](https://docs.apollo.io/):

- **People enrichment**(세 skill 모두 사용)은 person당 1 credit을 소비합니다
- **Bulk enrichment**(`/apollo:prospect`, `/apollo:sequence-load`)는 batch의 person당 1 credit을 소비합니다
- Plugin은 credit을 소비하기 전에 항상 경고합니다

## Credit

- **MCP Server** by [Apollo.io](https://docs.apollo.io/)
- **Plugin Specification** by [Anthropic](https://docs.anthropic.com/)

## 라이선스

MIT — 자세한 내용은 [LICENSE](LICENSE)를 참고하세요.
