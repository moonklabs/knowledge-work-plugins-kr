# Product Management 플러그인

Product management 플러그인입니다. 주로 Anthropic의 agentic 데스크톱 애플리케이션인 [Cowork](https://claude.com/product/cowork)용으로 설계되었지만 Claude Code에서도 동작합니다. Feature spec 작성, roadmap 관리, stakeholder communication, user research synthesis, competitor analysis, product metric tracking까지 PM workflow 전반을 다룹니다.

## 설치

```text
claude plugins add knowledge-work-plugins/product-management
```

## 주요 기능

이 플러그인은 다음을 도와주는 AI-powered product management partner를 제공합니다.

- **Feature Specs & PRDs** — Problem statement 또는 feature idea에서 structured PRD를 생성합니다. User story, requirement prioritization, success metric, scope management를 포함합니다.
- **Roadmap Planning** — Product roadmap을 create, update, reprioritize합니다. Now/Next/Later, quarterly theme, OKR-aligned format, dependency mapping을 지원합니다.
- **Stakeholder Updates** — Executive, engineering, customer 등 audience에 맞춘 status update를 생성합니다. Connected tool에서 context를 가져와 weekly update 부담을 줄입니다.
- **User Research Synthesis** — Interview note, survey data, support ticket을 structured insight로 바꿉니다. Theme을 식별하고 persona를 만들며 supporting evidence와 함께 opportunity area를 드러냅니다.
- **Competitive Analysis** — Competitor를 조사하고 feature comparison, positioning analysis, strategic implication이 포함된 brief를 생성합니다.
- **Metrics Review** — Product metric을 분석하고 trend를 식별하며 target과 비교해 actionable insight를 제공합니다.
- **Product Brainstorming** — Sharp sparring partner와 problem space를 탐색하고 idea를 생성하며 product thinking을 stress-test합니다. How Might We, Jobs-to-be-Done, First Principles, Opportunity Solution Trees 같은 framework를 활용합니다.

## 명령

| 명령 | 설명 |
|---|---|
| `/write-spec` | Problem statement에서 feature spec 또는 PRD를 작성합니다. |
| `/roadmap-update` | Roadmap을 update, create, reprioritize합니다. |
| `/stakeholder-update` | Weekly, monthly, launch용 stakeholder update를 생성합니다. |
| `/synthesize-research` | Interview, survey, ticket에서 user research를 종합합니다. |
| `/competitive-brief` | Competitive analysis brief를 생성합니다. |
| `/metrics-review` | Product metric을 review하고 analyze합니다. |
| `/brainstorm` | Product idea, problem space, strategic question을 sharp thinking partner와 brainstorm합니다. |

## 스킬

| 스킬 | 다루는 내용 |
|---|---|
| `feature-spec` | PRD structure, user story, requirement categorization, acceptance criteria를 다룹니다. |
| `roadmap-management` | RICE, MoSCoW 같은 prioritization framework, roadmap format, dependency mapping을 다룹니다. |
| `stakeholder-comms` | Audience별 update template, risk communication, decision documentation을 다룹니다. |
| `user-research-synthesis` | Thematic analysis, affinity mapping, persona development, opportunity sizing을 다룹니다. |
| `competitive-analysis` | Feature comparison matrix, positioning analysis, win/loss analysis를 다룹니다. |
| `metrics-tracking` | Product metric hierarchy, OKR goal setting, dashboard design, review cadence를 다룹니다. |
| `product-brainstorming` | Thinking partner로서 product idea를 brainstorm하고 problem space를 explore하며 assumption에 challenge합니다. New opportunity 탐색, product problem solution 생성, idea stress-test, PM이 방향을 정하기 전 sharp sparring partner와 생각을 정리해야 할 때 사용합니다. |

## 예시 워크플로

### PRD 작성

```text
You: /write-spec
Claude: 어떤 feature 또는 problem을 spec으로 만들까요?
You: Enterprise customer를 위해 SSO support를 추가해야 해
Claude: [Target user, constraint, success metric을 질문합니다]
Claude: [Problem statement, user story, requirement, success metric, open question이 포함된 full PRD를 생성합니다]
```

### Stakeholder update 준비

```text
You: /stakeholder-update
Claude: 어떤 update인가요? (weekly, monthly, launch, ad-hoc)
You: Exec team용 weekly update
Claude: [Project management tool, chat, docs에서 context를 가져옵니다]
Claude: [Progress, decision, risk, next milestone이 포함된 executive summary를 생성합니다]
```

### User research 종합

```text
You: /synthesize-research
Claude: 어떤 research를 종합할까요? Interview note를 paste하거나 file을 upload하거나 connected source에서 가져올 수 있습니다.
You: [8개 interview transcript를 붙여넣습니다]
Claude: [Interview 전반의 theme, pattern, insight를 식별합니다]
Claude: [Key finding, persona, opportunity area, evidence가 포함된 synthesis를 생성합니다]
```

### Product idea brainstorming

```text
You: /brainstorm 우리 product에 AI-powered search를 추가해야 할까?
Claude: AI search로 바로 가기 전에, 사용자가 지금 search에서 겪는 문제는 무엇인가요?
You: 원하는 걸 못 찾아. "X를 못 찾겠다"는 support ticket이 두 번째로 많아.
Claude: 이건 search technology 문제가 아니라 discovery problem일 수도 있습니다. 조금 더 파고들어볼게요...
Claude: [Issue가 search quality, information architecture, content discoverability 중 어디에 있는지 탐색합니다]
Claude: [Better taxonomy부터 AI search, guided navigation까지 5가지 접근을 생성합니다]
Claude: [Assumption에 challenge하고 가장 risky한 가정을 검증할 cheapest experiment를 제안합니다]
```

### Competitive analysis

```text
You: /competitive-brief
Claude: 어떤 competitor 또는 feature area를 분석할까요?
You: 우리 onboarding flow를 [competitor A], [competitor B]와 비교해줘
Claude: [Competitor onboarding approach를 조사합니다]
Claude: [Feature comparison, strength/weakness, strategic implication이 포함된 brief를 생성합니다]
```

## Data source

> 익숙하지 않은 placeholder가 보이거나 연결된 도구를 확인해야 한다면 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

가장 좋은 경험을 위해 project management 및 communication tool을 연결하세요. 연결하지 않았다면 context를 직접 제공하면 됩니다.

**포함된 MCP connection:**
- Chat(Slack): team context 및 stakeholder thread
- Project tracker(Linear, Asana, monday.com, ClickUp, Atlassian): roadmap integration, ticket context, status tracking
- Knowledge base(Notion): existing spec, research, meeting note
- Design(Figma): design context 및 handoff
- Product analytics(Amplitude, Pendo): usage data, metric, behavioral analysis
- User feedback(Intercom): support ticket, feature request, user conversation
- Meeting transcription(Fireflies): meeting note 및 discussion context

**추가 option:**
- Category별 alternative tool은 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.
