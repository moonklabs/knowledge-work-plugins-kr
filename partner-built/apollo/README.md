# Claude Code 및 Cowork용 Apollo 플러그인

[Apollo.io](https://www.apollo.io/)로 잠재 고객을 찾고, 리드를 보강하고, 아웃리치 시퀀스에 적재합니다. Apollo MCP Server 기반이며 **원클릭 통합**을 지원합니다.

---

## 원클릭 MCP server 통합

이 플러그인은 설치 시 **Apollo MCP Server를 자동으로 구성**합니다. 수동 서버 설정이나 config file 수정 없이 플러그인을 설치하고 Apollo 계정으로 인증하면 됩니다.

---

## 주요 스킬

이 플러그인은 여러 Apollo API를 완성된 작업 흐름으로 연결하는 고가치 스킬을 제공합니다.

| 스킬 | 수행 내용 |
|---|---|
| `/apollo:enrich-lead` | 이름, LinkedIn URL, email을 넣으면 email, phone, company intel, 다음 액션이 포함된 전체 contact card를 받습니다. |
| `/apollo:prospect` | ICP를 자연어로 설명하면 보강된 decision-maker lead의 순위표를 받습니다. |
| `/apollo:sequence-load` | 리드를 찾고 보강한 뒤 아웃리치 시퀀스에 일괄 적재합니다. 중복 제거와 등록을 처리합니다. |

### `/apollo:enrich-lead`

이름, 회사, LinkedIn URL, email을 넣으면 email, phone, 직함, 위치, 회사 세부정보, 추천 다음 액션이 포함된 완성된 contact card를 받습니다. "CEO of Figma" 같은 fuzzy lookup을 처리하고 exact match가 실패하면 search로 fallback합니다.

### `/apollo:prospect`

ICP를 자연어로 설명하세요. 파이프라인이 일치하는 회사를 검색하고 firmographic data를 일괄 보강하며 decision maker를 찾고 bulk enrichment로 contact info를 reveal한 뒤 ICP 적합도 점수가 포함된 리드 순위표를 반환합니다.

### `/apollo:sequence-load`

타기팅 기준에 맞는 contact를 찾고 보강한 뒤 중복 제거와 함께 contact로 생성하고 existing Apollo sequence에 bulk-add합니다. 등록 전에 후보를 미리 보고 이후 전체 요약을 보여줍니다.

---

## 설치

### Cowork

아래 링크를 클릭해 한 번에 설치합니다.

[Install in Cowork](https://claude.ai/desktop/customize/plugins/new?marketplace=apolloio/apollo-mcp-plugin&plugin=apollo)

MCP server가 올바르게 시작되도록 Cowork를 다시 시작하세요.

### Claude Code

#### 1. 이 plugin marketplace 추가

Claude Code에서 실행합니다.

```text
/plugin marketplace add apolloio/apollo-mcp-plugin
```

#### 2. Plugin 설치

```text
/plugin install apollo@apollo-plugin-marketplace
```

#### 3. Claude Code 재시작

이렇게 하면 MCP server가 올바르게 시작됩니다.

---

## 인증

Apollo MCP Server는 **OAuth**를 지원합니다.

1. 설치 후 Claude Code 또는 Cowork에서 `/mcp`를 실행합니다.
2. **Apollo** server를 선택하고 **Authenticate**를 클릭합니다.
3. Browser에서 Apollo.io 로그인을 완료합니다.
4. 완료되면 모든 명령을 사용할 준비가 됩니다.

---

## Apollo credit

일부 작업은 [Apollo credits](https://docs.apollo.io/)를 소비합니다.

- **People enrichment**(세 skill 모두 사용)는 사람당 1 credit을 소비합니다.
- **Bulk enrichment**(`/apollo:prospect`, `/apollo:sequence-load`)는 batch의 사람당 1 credit을 소비합니다.
- 플러그인은 credit을 소비하기 전에 항상 경고합니다.

## 크레딧

- **MCP Server** by [Apollo.io](https://docs.apollo.io/)
- **Plugin Specification** by [Anthropic](https://docs.anthropic.com/)

## 라이선스

MIT — 자세한 내용은 [LICENSE](LICENSE)를 참고하세요.
