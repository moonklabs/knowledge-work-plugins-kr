# Product Management 플러그인

Product Management 플러그인입니다. 주로 Anthropic의 agentic 데스크톱 애플리케이션인 [Cowork](https://claude.com/product/cowork)용으로 설계되었지만 Claude Code에서도 동작합니다. 기능 명세 작성, roadmap 관리, 이해관계자 커뮤니케이션, 사용자 리서치 종합, 경쟁사 분석, 제품 지표 추적까지 PM 워크플로 전반을 다룹니다.

## 설치

```text
claude plugins add knowledge-work-plugins/product-management
```

## 주요 기능

이 플러그인은 다음 작업을 돕는 AI 기반 Product Management 파트너를 제공합니다.

- **Feature Specs & PRDs** — Problem statement 또는 feature idea에서 구조화된 PRD를 생성합니다. 사용자 스토리, 요구사항 우선순위, 성공 지표, 범위 관리를 포함합니다.
- **Roadmap Planning** — Product roadmap을 생성, 갱신, 우선순위 재조정합니다. Now/Next/Later, 분기별 테마, OKR 정렬 형식, 의존성 매핑을 지원합니다.
- **Stakeholder Updates** — 임원진, 엔지니어링, 고객 등 대상에 맞춘 상태 업데이트를 생성합니다. 연결된 도구에서 맥락을 가져와 주간 업데이트 부담을 줄입니다.
- **User Research Synthesis** — 인터뷰 노트, 설문 데이터, 지원 티켓을 구조화된 인사이트로 바꿉니다. 주제를 식별하고 persona를 만들며 근거와 함께 기회 영역을 드러냅니다.
- **Competitive Analysis** — 경쟁사를 조사하고 기능 비교, 포지셔닝 분석, 전략적 시사점이 포함된 brief를 생성합니다.
- **Metrics Review** — 제품 지표를 분석하고 추세를 식별하며 목표와 비교해 실행 가능한 인사이트를 제공합니다.
- **Product Brainstorming** — 날카로운 sparring partner와 문제 공간을 탐색하고 아이디어를 생성하며 제품 사고를 stress-test합니다. How Might We, Jobs-to-be-Done, First Principles, Opportunity Solution Trees 같은 framework를 활용합니다.

## 명령

| 명령 | 설명 |
|---|---|
| `/write-spec` | Problem statement에서 기능 명세 또는 PRD를 작성합니다. |
| `/roadmap-update` | Roadmap을 갱신, 생성, 우선순위 재조정합니다. |
| `/stakeholder-update` | 주간, 월간, launch용 stakeholder update를 생성합니다. |
| `/synthesize-research` | 인터뷰, 설문, 티켓에서 사용자 리서치를 종합합니다. |
| `/competitive-brief` | 경쟁 분석 brief를 생성합니다. |
| `/metrics-review` | 제품 지표를 검토하고 분석합니다. |
| `/brainstorm` | 제품 아이디어, 문제 공간, 전략 질문을 날카로운 thinking partner와 brainstorm합니다. |

## 스킬

| 스킬 | 다루는 내용 |
|---|---|
| `feature-spec` | PRD 구조, 사용자 스토리, 요구사항 분류, 인수 기준을 다룹니다. |
| `roadmap-management` | RICE, MoSCoW 같은 우선순위 framework, roadmap 형식, 의존성 매핑을 다룹니다. |
| `stakeholder-comms` | 대상별 업데이트 템플릿, 위험 커뮤니케이션, 의사결정 문서화를 다룹니다. |
| `user-research-synthesis` | 주제 분석, affinity mapping, persona 개발, 기회 규모 산정을 다룹니다. |
| `competitive-analysis` | 기능 비교 매트릭스, 포지셔닝 분석, win/loss 분석을 다룹니다. |
| `metrics-tracking` | 제품 지표 계층, OKR 목표 설정, 대시보드 설계, 검토 주기를 다룹니다. |
| `product-brainstorming` | Thinking partner로서 제품 아이디어를 brainstorm하고 문제 공간을 explore하며 assumption에 challenge합니다. 새 기회 탐색, 제품 문제 해결안 생성, 아이디어 stress-test, PM이 방향을 정하기 전 날카로운 sparring partner와 생각을 정리해야 할 때 사용합니다. |

## 예시 워크플로

### PRD 작성

```text
You: /write-spec
Claude: 어떤 feature 또는 problem을 spec으로 만들까요?
You: Enterprise customer를 위해 SSO support를 추가해야 해
Claude: [대상 사용자, 제약 조건, 성공 지표를 질문합니다]
Claude: [Problem statement, 사용자 스토리, 요구사항, 성공 지표, open question이 포함된 전체 PRD를 생성합니다]
```

### Stakeholder update 준비

```text
You: /stakeholder-update
Claude: 어떤 update인가요? (weekly, monthly, launch, ad-hoc)
You: Exec team용 weekly update
Claude: [프로젝트 관리 도구, 채팅, 문서에서 맥락을 가져옵니다]
Claude: [진행 상황, 결정, 위험, 다음 milestone이 포함된 executive summary를 생성합니다]
```

### User research 종합

```text
You: /synthesize-research
Claude: 어떤 리서치를 종합할까요? 인터뷰 노트를 붙여넣거나 파일을 업로드하거나 연결된 출처에서 가져올 수 있습니다.
You: [8개 인터뷰 transcript를 붙여넣습니다]
Claude: [인터뷰 전반의 주제, 패턴, 인사이트를 식별합니다]
Claude: [핵심 발견, persona, 기회 영역, 근거가 포함된 synthesis를 생성합니다]
```

### Product idea brainstorming

```text
You: /brainstorm 우리 product에 AI 기반 search를 추가해야 할까?
Claude: AI search로 바로 가기 전에, 사용자가 지금 search에서 겪는 문제는 무엇인가요?
You: 원하는 걸 못 찾아. "X를 못 찾겠다"는 support ticket이 두 번째로 많아.
Claude: 이건 search technology 문제가 아니라 discovery problem일 수도 있습니다. 조금 더 파고들어볼게요...
Claude: [이슈가 검색 품질, 정보 구조, 콘텐츠 발견 가능성 중 어디에 있는지 탐색합니다]
Claude: [더 나은 taxonomy부터 AI search, guided navigation까지 5가지 접근을 생성합니다]
Claude: [가정에 challenge하고 가장 위험한 가정을 검증할 가장 저렴한 experiment를 제안합니다]
```

### Competitive analysis

```text
You: /competitive-brief
Claude: 어떤 competitor 또는 feature area를 분석할까요?
You: 우리 onboarding flow를 [competitor A], [competitor B]와 비교해줘
Claude: [경쟁사의 onboarding approach를 조사합니다]
Claude: [기능 비교, 강점/약점, 전략적 시사점이 포함된 brief를 생성합니다]
```

## 데이터 출처

> 익숙하지 않은 placeholder가 보이거나 연결된 도구를 확인해야 한다면 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

가장 좋은 경험을 위해 프로젝트 관리 및 커뮤니케이션 도구를 연결하세요. 연결하지 않았다면 맥락을 직접 제공하면 됩니다.

**포함된 MCP 연결:**
- 채팅(Slack): team context 및 stakeholder thread
- Project tracker(Linear, Asana, monday.com, ClickUp, Atlassian): roadmap integration, ticket context, status tracking
- Knowledge base(Notion): existing spec, research, meeting note
- Design(Figma): design context 및 handoff
- Product analytics(Amplitude, Pendo): usage data, metric, behavioral analysis
- User feedback(Intercom): support ticket, feature request, user conversation
- Meeting transcription(Fireflies): meeting note 및 discussion context

**추가 옵션:**
- 범주별 대체 도구는 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.
