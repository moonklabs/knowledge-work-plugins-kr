# Brand Voice 플러그인

Claude Cowork용 [Tribe AI](https://tribe.ai) plugin입니다. Cowork launch partner로 제작되었습니다.

회사를 알아볼 수 있게 만드는 brand knowledge는 보통 편리한 한곳에 있지 않습니다. 2022년 deck, 마지막 rebrand 이후 아무도 업데이트하지 않은 Confluence page, 오래 일한 senior 몇 명의 감각 속에 흩어져 있습니다. Sales rep이 AI로 outreach를 만들고 new hire가 첫 주부터 content를 만들 때, 바로 그 지식이 사라집니다.

Brand Voice는 흩어진 brand material을 enforce 가능한 AI guardrail로 바꿉니다. Confluence, Google Drive, Box, SharePoint, Slack, Gong, Granola를 검색해 회사가 실제로 어떻게 소통하는지 찾아낸 뒤, LLM-ready brand guideline을 만들고 AI-generated content를 그 기준으로 검증합니다. Claude는 더 빠르게 쓰는 데서 그치지 않고, 여러분답게 씁니다.

## 기능

### 1. Brand discovery
Brand knowledge는 Notion, Confluence, Google Drive, Gong, Slack, 수년간의 sales call 및 meeting transcript에 흩어져 있습니다. Brand Voice는 style guide, pitch deck, email template, transcript, design system을 모두 검색해 가장 강한 brand signal을 하나의 최신 source of truth로 정리합니다. 3년 전 style guide의 문구가 아니라, 실제로 가장 잘하는 사람들이 어떻게 소통하는지에 기반합니다.

**Slash 명령:** `/brand-voice:discover-brand`

```
/brand-voice:discover-brand
/brand-voice:discover-brand Acme Corp
```

### 2. Guideline generation
자료를 LLM-ready guideline으로 종합합니다. Voice pillar, tone parameter, Claude에게 명확한 operating boundary를 주는 "We Are / We Are Not" framework, 실제 company language를 반영하는 terminology standard를 만듭니다. Veteran team을 on-brand로 유지하는 같은 guardrail 덕분에 new hire도 세 달째가 아니라 첫 주부터 quality content를 만들 수 있습니다.

**Slash 명령:** `/brand-voice:generate-guidelines`

```
/brand-voice:generate-guidelines
/brand-voice:generate-guidelines from the discovery report and these 3 PDFs
```

### 3. Brand voice enforcement
Sales email, proposal, marketing page, press release 같은 모든 AI-generated content는 처음부터 guideline 기준으로 작성됩니다. Voice는 일정하게 유지하고 tone은 context에 맞게 유연하게 바뀝니다. Cold email, enterprise proposal, LinkedIn post에 따라 formality, energy, technical depth가 자동 조정됩니다. Tone drift와 positioning gap은 prospect나 investor에게 도달하기 전에 잡아냅니다.

**Slash 명령:** `/brand-voice:enforce-voice`

```
/brand-voice:enforce-voice Draft a cold email to a VP of Sales at a mid-market SaaS company
/brand-voice:enforce-voice Write a LinkedIn post announcing our new feature
```

### Open question
Plugin이 conflicting document, missing guideline, stated brand와 practiced brand의 차이처럼 해결할 수 없는 ambiguity를 만나면 team discussion용 open question을 표시합니다. 모든 question에는 agent recommendation이 포함되어 ambiguity를 dead end가 아니라 "confirm or override" interaction으로 바꿉니다.

## MCP connector

| Connector | URL | 목적 |
|-----------|-----|---------|
| **Notion** | `https://mcp.notion.com/mcp` | Discovery backbone — federates across connected Google Drive, SharePoint, OneDrive, Slack, Jira. Also stores output guidelines. |
| **Atlassian** | `https://mcp.atlassian.com/v1/mcp` | Deep Confluence search + Jira context for Atlassian-heavy enterprises |
| **Box** | `https://mcp.box.com` | Cloud file storage — official brand docs, shared decks, and style guides often live here |
| **Microsoft 365** | `https://microsoft365.mcp.claude.com/mcp` | SharePoint, OneDrive, Outlook, Teams — enterprise document storage and email templates |
| **Figma** | `https://mcp.figma.com/mcp` | Brand design systems — color, typography, design tokens inform voice |
| **Gong** | `https://mcp.gong.io/mcp` | Enterprise conversation intelligence — sales call transcripts and analysis |
| **Granola** | `https://mcp.granola.ai/mcp` | Meeting intelligence — transcripts and notes from sales, customer, and strategy meetings |

### Native integration

이 platform들은 native Claude integration입니다. MCP connector install이 필요 없습니다. 사용자가 Claude Desktop 또는 Cowork에서 연결하면 tool로 사용할 수 있습니다.

| Integration | 목적 |
|-------------|---------|
| **Google Drive** | Shared brand documents, style guides, marketing materials, Google Docs and Slides |
| **Slack** | Brand discussions, channel searches, pinned brand guidelines, informal voice patterns |

## 빠른 시작

1. Plugin을 설치하고 Claude Cowork를 엽니다.
2. 최소 하나의 platform을 연결합니다. Notion을 권장합니다. Google Drive, SharePoint, Slack, Jira를 federate합니다.
3. `/brand-voice:discover-brand`를 실행합니다. Claude가 connected knowledge base에서 brand material을 자동 검색합니다.
4. `/brand-voice:generate-guidelines`를 실행해 discovery report 기반의 durable guideline을 만듭니다.
5. Sales email, proposal, LinkedIn post 등 customer-facing content를 만들 때 `/brand-voice:enforce-voice`를 사용합니다.

원한다면 Claude에게 specific document를 직접 지정할 수도 있습니다. 어느 쪽이든 plugin이 process를 안내합니다.

Brand Voice는 현재 individual level에서 동작하며, team-wide enforcement는 곧 제공될 예정입니다.

### Project별 설정

`settings/brand-voice.local.md.example`을 project의 `.claude/brand-voice.local.md`로 복사하고 company name, enabled platform, known brand material location을 채우세요.

## 파일 구조

```
├── .claude-plugin/
│   └── plugin.json                              # Plugin manifest
├── .mcp.json                                    # 7 MCP server connections
├── README.md
├── agents/
│   ├── discover-brand.md                        # Autonomous platform search agent
│   ├── content-generation.md                    # Brand-aligned content creation
│   ├── conversation-analysis.md                 # Sales call transcript analysis
│   ├── document-analysis.md                     # Brand document parsing
│   └── quality-assurance.md                     # Compliance and open questions audit
├── commands/
│   ├── discover-brand.md                        # /brand-voice:discover-brand
│   ├── enforce-voice.md                         # /brand-voice:enforce-voice
│   └── generate-guidelines.md                   # /brand-voice:generate-guidelines
├── settings/
│   └── brand-voice.local.md.example             # Per-project settings template
└── skills/
    ├── discover-brand/
    │   ├── SKILL.md                             # Discovery orchestration
    │   └── references/
    │       ├── search-strategies.md             # Platform-specific query patterns
    │       └── source-ranking.md                # Ranking algorithm and categories
    ├── brand-voice-enforcement/
    │   ├── SKILL.md                             # Enforcement orchestration
    │   └── references/
    │       ├── before-after-examples.md         # Content type transformation examples
    │       └── voice-constant-tone-flexes.md    # "We Are / We Are Not" + tone matrix
    └── guideline-generation/
        ├── SKILL.md                             # Generation orchestration
        └── references/
            ├── confidence-scoring.md            # Scoring methodology
            └── guideline-template.md            # Full output template
```

## 아키텍처

**Skills**는 domain knowledge를 제공하고 workflow를 orchestrate합니다. User intent에 따라 자동 활성화됩니다.

**Agents**는 platform search, document analysis, transcript parsing, content generation, quality validation 같은 heavy autonomous work를 처리합니다.

**Commands**는 skill workflow를 trigger하는 explicit user entry point입니다.

**핵심 설계 결정:**
- Voice는 일정하고 tone은 유연합니다. Enforcement를 위한 명확한 mental model입니다.
- Discovery agent는 autonomous하지만 accountable합니다. Provenance와 conflict를 함께 보여줍니다.
- Open question은 failure가 아니라 feature입니다. 모든 ambiguity에는 recommendation이 포함됩니다.
- Progressive disclosure를 사용합니다. Frontmatter는 lean하게, `SKILL.md`는 focused하게 유지하고 detail은 `references/`에 둡니다.
- Notion AI Search를 federated discovery engine으로 사용합니다. 하나의 API가 connected source를 통해 8개 이상 platform을 검색합니다.
- Google Drive와 Slack은 native Claude integration입니다. MCP connector가 필요 없습니다.
