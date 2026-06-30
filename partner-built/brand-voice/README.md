# Brand Voice 플러그인

Claude Cowork용 [Tribe AI](https://tribe.ai) 플러그인입니다. Cowork 출시 파트너로 제작되었습니다.

회사를 알아볼 수 있게 만드는 브랜드 지식은 보통 편리한 한곳에 있지 않습니다. 2022년 덱, 마지막 리브랜딩 이후 아무도 업데이트하지 않은 Confluence 페이지, 오래 일한 선임 몇 명의 감각 속에 흩어져 있습니다. 영업 담당자가 AI로 아웃리치를 만들고 신규 입사자가 첫 주부터 콘텐츠를 만들 때, 바로 그 지식이 사라집니다.

Brand Voice는 흩어진 브랜드 자료를 적용 가능한 AI 가드레일로 바꿉니다. Confluence, Google Drive, Box, SharePoint, Slack, Gong, Granola를 검색해 회사가 실제로 어떻게 소통하는지 찾아낸 뒤, LLM에 바로 쓸 수 있는 브랜드 가이드라인을 만들고 AI 생성 콘텐츠를 그 기준으로 검증합니다. Claude는 더 빠르게 쓰는 데서 그치지 않고, 여러분답게 씁니다.

## 기능

### 1. 브랜드 발견
브랜드 지식은 Notion, Confluence, Google Drive, Gong, Slack, 수년간의 영업 통화와 미팅 전사에 흩어져 있습니다. Brand Voice는 스타일 가이드, 피치 덱, 이메일 템플릿, 전사, 디자인 시스템을 모두 검색해 가장 강한 브랜드 신호를 하나의 최신 기준점으로 정리합니다. 3년 전 스타일 가이드의 문구가 아니라, 실제로 가장 잘하는 사람들이 어떻게 소통하는지에 기반합니다.

**슬래시 명령:** `/brand-voice:discover-brand`

```
/brand-voice:discover-brand
/brand-voice:discover-brand Acme Corp
```

### 2. 가이드라인 생성
자료를 LLM에 바로 쓸 수 있는 가이드라인으로 종합합니다. 보이스 축, 톤 매개변수, Claude에게 명확한 운영 경계를 주는 "We Are / We Are Not" framework, 실제 회사 언어를 반영하는 용어 기준을 만듭니다. 숙련된 팀을 브랜드에 맞게 유지하는 같은 가드레일 덕분에 신입도 세 달째가 아니라 첫 주부터 품질 높은 콘텐츠를 만들 수 있습니다.

**슬래시 명령:** `/brand-voice:generate-guidelines`

```
/brand-voice:generate-guidelines
/brand-voice:generate-guidelines discovery report와 이 PDF 3개를 바탕으로 만들어줘
```

### 3. Brand voice 적용
영업 이메일, 제안서, 마케팅 페이지, 보도자료 같은 모든 AI 생성 콘텐츠는 처음부터 가이드라인 기준으로 작성됩니다. Voice는 일정하게 유지하고 tone은 맥락에 맞게 유연하게 바뀝니다. Cold email, enterprise proposal, LinkedIn post에 따라 격식, 에너지, 기술적 깊이가 자동 조정됩니다. 톤 이탈과 포지셔닝 공백은 잠재 고객이나 투자자에게 도달하기 전에 잡아냅니다.

**슬래시 명령:** `/brand-voice:enforce-voice`

```
/brand-voice:enforce-voice mid-market SaaS 회사의 VP of Sales에게 보낼 cold email을 작성해줘
/brand-voice:enforce-voice 새 기능 출시를 알리는 LinkedIn post를 작성해줘
```

### 미해결 질문
플러그인이 서로 충돌하는 문서, 누락된 가이드라인, 선언된 브랜드와 실제 사용되는 브랜드의 차이처럼 해결할 수 없는 모호성을 만나면 팀 논의용 미해결 질문을 표시합니다. 모든 질문에는 에이전트 추천이 포함되어 모호성을 막다른 길이 아니라 "confirm or override" 상호작용으로 바꿉니다.

## MCP 커넥터

| 커넥터 | URL | 목적 |
|-----------|-----|---------|
| **Notion** | `https://mcp.notion.com/mcp` | Discovery backbone — 연결된 Google Drive, SharePoint, OneDrive, Slack, Jira를 통합 검색합니다. 생성된 가이드라인도 저장합니다. |
| **Atlassian** | `https://mcp.atlassian.com/v1/mcp` | Atlassian 중심 조직을 위한 깊은 Confluence 검색 + Jira 맥락 |
| **Box** | `https://mcp.box.com` | 클라우드 파일 스토리지 — 공식 브랜드 문서, 공유 덱, style guide가 자주 있는 곳 |
| **Microsoft 365** | `https://microsoft365.mcp.claude.com/mcp` | SharePoint, OneDrive, Outlook, Teams — 엔터프라이즈 문서 저장소와 이메일 템플릿 |
| **Figma** | `https://mcp.figma.com/mcp` | 브랜드 디자인 시스템 — 색상, 타이포그래피, 디자인 토큰이 voice 판단에 도움을 줍니다 |
| **Gong** | `https://mcp.gong.io/mcp` | 엔터프라이즈 대화 인텔리전스 — 영업 통화 전사와 분석 |
| **Granola** | `https://mcp.granola.ai/mcp` | 미팅 인텔리전스 — 영업, 고객, 전략 미팅의 전사와 노트 |

### 네이티브 통합

이 플랫폼들은 Claude 기본 통합입니다. MCP 커넥터 설치가 필요 없습니다. 사용자가 Claude Desktop 또는 Cowork에서 연결하면 도구로 사용할 수 있습니다.

| 통합 | 목적 |
|-------------|---------|
| **Google Drive** | 공유 브랜드 문서, style guide, 마케팅 자료, Google Docs 및 Slides |
| **Slack** | 브랜드 논의, 채널 검색, 고정된 브랜드 가이드라인, 비공식 voice pattern |

## 빠른 시작

1. 플러그인을 설치하고 Claude Cowork를 엽니다.
2. 최소 하나의 플랫폼을 연결합니다. Notion을 권장합니다. Google Drive, SharePoint, Slack, Jira를 통합 검색합니다.
3. `/brand-voice:discover-brand`를 실행합니다. Claude가 연결된 지식 베이스에서 브랜드 자료를 자동 검색합니다.
4. `/brand-voice:generate-guidelines`를 실행해 discovery report 기반의 지속 가능한 guideline을 만듭니다.
5. 영업 이메일, 제안서, LinkedIn post 등 고객 대상 콘텐츠를 만들 때 `/brand-voice:enforce-voice`를 사용합니다.

원한다면 Claude에게 특정 문서를 직접 지정할 수도 있습니다. 어느 쪽이든 플러그인이 과정을 안내합니다.

Brand Voice는 현재 개인 단위에서 동작하며, 팀 전체 적용은 곧 제공될 예정입니다.

### 프로젝트별 설정

`settings/brand-voice.local.md.example`을 project의 `.claude/brand-voice.local.md`로 복사하고 회사명, 활성화할 플랫폼, 알려진 브랜드 자료 위치를 채우세요.

## 파일 구조

```
├── .claude-plugin/
│   └── plugin.json                              # 플러그인 매니페스트
├── .mcp.json                                    # 7개 MCP 서버 연결
├── README.md
├── agents/
│   ├── discover-brand.md                        # 자율 platform 검색 agent
│   ├── content-generation.md                    # 브랜드 정렬 콘텐츠 생성
│   ├── conversation-analysis.md                 # 영업 통화 전사 분석
│   ├── document-analysis.md                     # 브랜드 문서 파싱
│   └── quality-assurance.md                     # 컴플라이언스 및 open question audit
├── commands/
│   ├── discover-brand.md                        # /brand-voice:discover-brand
│   ├── enforce-voice.md                         # /brand-voice:enforce-voice
│   └── generate-guidelines.md                   # /brand-voice:generate-guidelines
├── settings/
│   └── brand-voice.local.md.example             # 프로젝트별 설정 template
└── skills/
    ├── discover-brand/
    │   ├── SKILL.md                             # Discovery 흐름 조율
    │   └── references/
    │       ├── search-strategies.md             # Platform별 query pattern
    │       └── source-ranking.md                # Ranking algorithm 및 category
    ├── brand-voice-enforcement/
    │   ├── SKILL.md                             # Enforcement 흐름 조율
    │   └── references/
    │       ├── before-after-examples.md         # Content type별 transformation 예시
    │       └── voice-constant-tone-flexes.md    # "We Are / We Are Not" + tone matrix
    └── guideline-generation/
        ├── SKILL.md                             # Generation 흐름 조율
        └── references/
            ├── confidence-scoring.md            # Scoring 방법론
            └── guideline-template.md            # 전체 output template
```

## 아키텍처

**스킬**은 도메인 지식과 워크플로 조율을 제공합니다. 사용자 의도에 따라 자동 활성화됩니다.

**에이전트**는 플랫폼 검색, 문서 분석, 전사 파싱, 콘텐츠 생성, 품질 검증 같은 무거운 자율 작업을 처리합니다.

**명령**은 스킬 워크플로를 트리거하는 명시적 사용자 진입점입니다.

**핵심 설계 결정:**
- Voice는 일정하고 tone은 유연합니다. Enforcement를 위한 명확한 mental model입니다.
- Discovery agent는 자율적으로 동작하지만 책임 추적이 가능합니다. 출처와 충돌을 함께 보여줍니다.
- Open question은 실패가 아니라 기능입니다. 모든 모호성에는 추천이 포함됩니다.
- 점진적 공개를 사용합니다. Frontmatter는 간결하게, `SKILL.md`는 집중된 상태로 유지하고 세부사항은 `references/`에 둡니다.
- Notion AI Search를 연합 discovery engine으로 사용합니다. 하나의 API가 연결된 source를 통해 8개 이상 플랫폼을 검색합니다.
- Google Drive와 Slack은 native Claude integration입니다. MCP 커넥터가 필요 없습니다.
