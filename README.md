# 지식 업무 플러그인 (한글버전)

Claude를 사용자의 역할, 팀, 회사에 맞는 전문가로 바꿔 주는 플러그인 모음입니다. [Claude Cowork](https://claude.com/product/cowork)용으로 제작되었으며 [Claude Code](https://claude.com/product/claude-code)와도 호환됩니다.

## 플러그인이 필요한 이유

Cowork에서는 목표를 정하면 Claude가 완성도 높은 전문 작업물을 만들어 줍니다. 플러그인을 사용하면 한 단계 더 나아가, Claude에게 선호하는 작업 방식, 가져와야 할 도구와 데이터, 중요한 워크플로우 처리 방식, 노출할 슬래시 명령을 알려 줄 수 있습니다. 그래서 팀은 더 일관되고 품질 높은 결과를 얻습니다.

각 플러그인은 특정 직무에 필요한 스킬, 커넥터, 슬래시 명령, 하위 에이전트를 묶습니다. 기본 상태에서도 해당 역할을 돕기 위한 강한 출발점을 제공하지만, 진짜 가치는 회사의 도구, 용어, 프로세스에 맞게 커스터마이즈할 때 나옵니다. 그렇게 하면 Claude가 처음부터 팀을 위해 만들어진 것처럼 일합니다.

## 플러그인 마켓플레이스

우리는 내부 업무에서 만들고 영감을 얻은 11개의 플러그인을 오픈소스로 제공합니다.

| 플러그인 | 도움을 주는 일 | 커넥터 |
|--------|-------------|------------|
| **[productivity](./productivity)** | 작업, 캘린더, 일일 워크플로, 개인 컨텍스트를 관리해 반복 설명 시간을 줄입니다. | Slack, Notion, Asana, Linear, Jira, Monday, ClickUp, Microsoft 365 |
| **[sales](./sales)** | 잠재 고객 조사, 통화 준비, 파이프라인 검토, 아웃리치 작성, 경쟁 배틀카드 작성을 돕습니다. | Slack, HubSpot, Close, Clay, ZoomInfo, Notion, Jira, Fireflies, Microsoft 365 |
| **[customer-support](./customer-support)** | 티켓 분류, 답변 초안, 에스컬레이션 패키징, 고객 컨텍스트 조사, 해결 이슈의 지식 베이스화를 지원합니다. | Slack, Intercom, HubSpot, Guru, Jira, Notion, Microsoft 365 |
| **[product-management](./product-management)** | 스펙 작성, 로드맵 계획, 사용자 리서치 종합, 이해관계자 업데이트, 경쟁 환경 추적을 돕습니다. | Slack, Linear, Asana, Monday, ClickUp, Jira, Notion, Figma, Amplitude, Pendo, Intercom, Fireflies |
| **[marketing](./marketing)** | 콘텐츠 초안, 캠페인 계획, 브랜드 보이스 적용, 경쟁사 브리핑, 채널별 성과 리포팅을 지원합니다. | Slack, Canva, Figma, HubSpot, Amplitude, Notion, Ahrefs, SimilarWeb, Klaviyo |
| **[legal](./legal)** | 계약 검토, NDA 분류, 컴플라이언스 검토, 위험 평가, 회의 준비, 템플릿 응답 작성을 돕습니다. | Slack, Box, Egnyte, Jira, Microsoft 365 |
| **[finance](./finance)** | 분개 준비, 계정 조정, 재무제표 생성, 차이 분석, 결산 관리, 감사 지원을 돕습니다. | Snowflake, Databricks, BigQuery, Slack, Microsoft 365 |
| **[data](./data)** | 데이터셋 질의, 시각화, 해석을 지원합니다. SQL 작성, 통계 분석, 대시보드 구축, 공유 전 검증을 돕습니다. | Snowflake, Databricks, BigQuery, Definite, Hex, Amplitude, Jira |
| **[enterprise-search](./enterprise-search)** | 이메일, 채팅, 문서, 위키 전반에서 필요한 정보를 찾습니다. 한 번의 질의로 회사 도구 전체를 검색합니다. | Slack, Notion, Guru, Jira, Asana, Microsoft 365 |
| **[bio-research](./bio-research)** | 문헌 검색, 유전체 분석, 표적 우선순위화 등 전임상 연구 도구와 데이터베이스를 연결해 초기 생명과학 R&D를 가속합니다. | PubMed, BioRender, bioRxiv, ClinicalTrials.gov, ChEMBL, Synapse, Wiley, Owkin, Open Targets, Benchling |
| **[cowork-plugin-management](./cowork-plugin-management)** | 새 플러그인을 만들거나 조직의 도구와 워크플로에 맞게 기존 플러그인을 커스터마이즈합니다. | - |

이 플러그인들은 Cowork에서 바로 설치하거나, GitHub에서 전체 컬렉션을 살펴보거나, 직접 새 플러그인을 만드는 출발점으로 사용할 수 있습니다.

## 현지화 릴리스 노트

### 2026-07-09

- upstream marketplace 업데이트를 `f96c57c`까지 동기화했습니다.
- 새 Cowork marketplace 항목 `Honeycomb`, `B12`, `SigNoz`의 설명을 한국어로 번역했습니다.
- 원격 `korean` 브랜치의 README 업데이트를 보존해 병합했습니다.
- marketplace 설명의 영어-only 잔여 항목과 `English (한국어)` 표시명 규칙을 검증했습니다.

## 시작하기

### Cowork

[claude.com/plugins](https://claude.com/plugins/)에서 플러그인을 설치합니다.

### Claude Code

```bash
# 먼저 마켓플레이스를 추가합니다
claude plugin marketplace add moonklabs/knowledge-work-plugins-kr

# 그다음 필요한 플러그인을 설치합니다
claude plugin install sales@knowledge-work-plugins-kr
```

설치가 끝나면 플러그인은 자동으로 활성화됩니다. 스킬은 관련 상황에서 자동으로 작동하고, 슬래시 명령은 세션에서 사용할 수 있습니다. 예: `/sales:call-prep`, `/data:write-query`.

## 플러그인 동작 방식

모든 플러그인은 같은 구조를 따릅니다.

```
plugin-name/
├── .claude-plugin/plugin.json   # 매니페스트
├── .mcp.json                    # 도구 연결
├── commands/                    # 명시적으로 실행하는 슬래시 명령
└── skills/                      # Claude가 자동으로 활용하는 도메인 지식
```

- **스킬**은 Claude가 유용한 도움을 주는 데 필요한 도메인 전문성, 모범 사례, 단계별 워크플로를 담습니다. Claude는 관련 상황에서 이를 자동으로 활용합니다.
- **명령**은 사용자가 명시적으로 실행하는 동작입니다. 예: `/finance:reconciliation`, `/product-management:write-spec`.
- **커넥터**는 Claude를 CRM, 프로젝트 추적기, 데이터 웨어하우스, 디자인 도구 같은 외부 도구에 연결합니다. 연결은 [MCP 서버](https://modelcontextprotocol.io/)를 통해 이뤄집니다.

모든 구성요소는 파일 기반입니다. Markdown과 JSON만 있으면 되며, 별도의 코드, 인프라, 빌드 단계가 없습니다.

## 우리 방식에 맞게 만들기

이 플러그인들은 범용 출발점입니다. 회사의 실제 업무 방식에 맞춰 조정할수록 훨씬 더 유용해집니다.

- **커넥터 교체** — `.mcp.json`을 수정해 회사의 실제 도구 스택을 가리키게 합니다.
- **회사 컨텍스트 추가** — 용어, 조직 구조, 프로세스를 스킬 파일에 넣어 Claude가 우리 업무 방식을 이해하게 합니다.
- **워크플로우 조정** — 교과서식 절차가 아니라 실제 팀이 일하는 방식에 맞게 스킬 지침을 수정합니다.
- **새 플러그인 만들기** — `cowork-plugin-management` 플러그인을 사용하거나 위 구조를 따라 아직 다루지 않은 역할과 워크플로우용 플러그인을 만듭니다.

팀이 플러그인을 만들고 공유할수록 Claude는 더 넓은 직무를 이해하는 다기능 전문가가 됩니다. 정의한 컨텍스트는 관련 상호작용마다 반영되므로, 리더와 관리자는 프로세스를 반복 설명하는 대신 프로세스를 개선하는 데 더 많은 시간을 쓸 수 있습니다.

## 기여하기

플러그인은 Markdown 파일입니다. 저장소를 fork하고, 변경한 뒤 PR을 제출하세요.
