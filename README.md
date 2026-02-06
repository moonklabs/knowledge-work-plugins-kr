# 지식 업무 플러그인 by Anthropic

Claude를 당신의 역할, 팀, 회사를 위한 전문가로 만들어주는 플러그인입니다. [Claude Cowork](https://claude.com/product/cowork)를 위해 만들어졌으며, [Claude Code](https://claude.com/product/claude-code)와도 호환됩니다.

## 플러그인이 필요한 이유

Cowork에서는 목표를 설정하면 Claude가 완성도 높은 결과물을 제공합니다. 플러그인을 사용하면 더 나아갈 수 있습니다. 원하는 작업 방식, 불러올 도구와 데이터, 중요한 워크플로 처리 방식, 노출할 슬래시 커맨드를 지정해 팀의 결과를 더 일관되고 더 좋게 만들 수 있습니다.

각 플러그인은 특정 직무에 필요한 스킬, 커넥터, 슬래시 커맨드, 하위 에이전트를 묶어서 제공합니다. 기본 제공 상태만으로도 해당 역할을 돕는 데 충분한 출발점을 제공합니다. 진짜 힘은 회사에 맞게 커스터마이즈할 때 나옵니다. 사용 중인 도구, 용어, 프로세스를 반영하면 Claude가 팀을 위해 만들어진 것처럼 동작합니다.

## 플러그인 마켓플레이스

우리의 작업에서 출발해 제작한 11개의 플러그인을 오픈 소스로 제공합니다.

| 플러그인 | 어떤 도움을 주는가 | 커넥터 |
|--------|-------------|------------|
| **[productivity](./productivity)** | 반복 업무 시간을 줄이기 위해 할 일, 일정, 일상 워크플로, 개인 컨텍스트를 관리합니다. | Slack, Notion, Asana, Linear, Jira, Monday, ClickUp, Microsoft 365 |
| **[sales](./sales)** | 잠재 고객 조사, 콜 준비, 파이프라인 리뷰, 아웃리치 초안 작성, 경쟁사 배틀카드 제작을 지원합니다. | Slack, HubSpot, Close, Clay, ZoomInfo, Notion, Jira, Fireflies, Microsoft 365 |
| **[customer-support](./customer-support)** | 티켓 분류, 응답 초안 작성, 에스컬레이션 패키징, 고객 컨텍스트 조사, 해결 이슈의 지식베이스화까지 돕습니다. | Slack, Intercom, HubSpot, Guru, Jira, Notion, Microsoft 365 |
| **[product-management](./product-management)** | 스펙 작성, 로드맵 계획, 사용자 리서치 합성, 이해관계자 업데이트, 경쟁 환경 추적을 지원합니다. | Slack, Linear, Asana, Monday, ClickUp, Jira, Notion, Figma, Amplitude, Pendo, Intercom, Fireflies |
| **[marketing](./marketing)** | 콘텐츠 초안 작성, 캠페인 계획, 브랜드 보이스 유지, 경쟁사 브리핑, 채널 성과 보고를 돕습니다. | Slack, Canva, Figma, HubSpot, Amplitude, Notion, Ahrefs, SimilarWeb, Klaviyo |
| **[legal](./legal)** | 계약 검토, NDA 분류, 컴플라이언스 대응, 리스크 평가, 미팅 준비, 템플릿 응답 작성에 도움을 줍니다. | Slack, Box, Egnyte, Jira, Microsoft 365 |
| **[finance](./finance)** | 전표 준비, 계정 조정, 재무제표 생성, 변동 분석, 마감 관리, 감사 지원을 돕습니다. | Snowflake, Databricks, BigQuery, Slack, Microsoft 365 |
| **[data](./data)** | SQL 작성, 통계 분석, 대시보드 구축, 공유 전 검증까지 데이터 작업 전반을 지원합니다. | Snowflake, Databricks, BigQuery, Hex, Amplitude, Jira |
| **[enterprise-search](./enterprise-search)** | 이메일, 채팅, 문서, 위키 등 회사 도구 전반을 하나의 쿼리로 검색합니다. | Slack, Notion, Guru, Jira, Asana, Microsoft 365 |
| **[bio-research](./bio-research)** | 문헌 검색, 유전체 분석, 타깃 우선순위화 등 전임상 연구 도구와 데이터베이스를 연결해 생명과학 R&D를 가속합니다. | PubMed, BioRender, bioRxiv, ClinicalTrials.gov, ChEMBL, Synapse, Wiley, Owkin, Open Targets, Benchling |
| **[cowork-plugin-management](./cowork-plugin-management)** | 조직별 도구와 워크플로에 맞게 새 플러그인을 만들거나 기존 플러그인을 커스터마이즈합니다. | — |

Cowork에서 바로 설치하거나, GitHub에서 전체 컬렉션을 확인하거나, 직접 만들어 사용할 수 있습니다.

## 시작하기

### Cowork

[claude.com/plugins](https://claude.com/plugins/)에서 플러그인을 설치하세요.

### Claude Code

```bash
# 먼저 마켓플레이스를 추가합니다.
claude plugin marketplace add anthropics/knowledge-work-plugins

# 그런 다음 특정 플러그인을 설치합니다.
claude plugin install sales@knowledge-work-plugins
```

설치 후 플러그인은 자동으로 활성화됩니다. 관련 상황에서 스킬이 자동으로 동작하며, 세션에서 슬래시 커맨드를 사용할 수 있습니다. 예: `/sales:call-prep`, `/data:write-query`.

## 플러그인 동작 방식

모든 플러그인은 동일한 구조를 따릅니다.

```
plugin-name/
├── .claude-plugin/plugin.json   # 매니페스트
├── .mcp.json                    # 도구 연결
├── commands/                    # 명시적으로 호출하는 슬래시 커맨드
└── skills/                      # Claude가 자동으로 사용하는 도메인 지식
```

- **스킬**은 Claude가 유용한 도움을 제공하기 위해 필요한 도메인 전문성, 모범 사례, 단계별 워크플로를 담습니다. 관련 상황에서 자동으로 사용됩니다.
- **커맨드**는 사용자가 명시적으로 트리거하는 작업입니다. 예: `/finance:reconciliation`, `/product-management:write-spec`.
- **커넥터**는 MCP 서버를 통해 CRM, 프로젝트 트래커, 데이터 웨어하우스, 디자인 도구 등 업무에 필요한 외부 도구와 Claude를 연결합니다.

모든 구성 요소는 파일 기반입니다. Markdown과 JSON만 사용하며, 코드나 인프라, 빌드 단계가 필요하지 않습니다.

## 우리 조직에 맞게 만들기

이 플러그인들은 범용 시작점입니다. 회사 실제 업무 방식에 맞게 커스터마이즈하면 훨씬 유용해집니다.

- **커넥터 교체**: `.mcp.json`을 편집해 조직의 도구 스택을 연결합니다.
- **회사 컨텍스트 추가**: 용어, 조직 구조, 프로세스를 스킬 파일에 추가해 Claude가 조직을 이해하도록 합니다.
- **워크플로 조정**: 실제 업무 방식에 맞게 스킬 지침을 조정합니다.
- **새 플러그인 만들기**: `cowork-plugin-management` 플러그인을 사용하거나 위 구조를 따라 새로운 역할용 플러그인을 만듭니다.

팀이 플러그인을 만들고 공유할수록 Claude는 크로스 펑셔널 전문가가 됩니다. 정의한 컨텍스트가 모든 관련 상호작용에 반영되어, 리더와 운영자는 프로세스 강제보다 개선에 더 많은 시간을 쓸 수 있습니다.

## 기여하기

플러그인은 Markdown 파일로 구성됩니다. 저장소를 포크하고 변경한 뒤 PR을 보내 주세요.
