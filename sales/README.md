# Sales 플러그인

Sales productivity plugin입니다. 주로 Anthropic의 agentic desktop application인 [Cowork](https://claude.com/product/cowork)용으로 설계되었지만 Claude Code에서도 동작합니다. Prospecting, outreach, pipeline management, call preparation, deal strategy를 돕습니다. 어떤 sales team에서도 사용할 수 있으며, web search와 사용자 입력만으로 standalone 동작하고 CRM, email 및 다른 tool을 연결하면 더 강력해집니다.

## 설치

```bash
claude plugins add knowledge-work-plugins/sales
```

## 명령

Slash command로 호출하는 명시적 workflow입니다:

| 명령 | 설명 |
|---|---|
| `/call-summary` | Process call notes or transcript — extract action items, draft follow-up, generate internal summary |
| `/forecast` | Generate a weighted sales forecast — upload CSV or describe your pipeline, set quota, get projections |
| `/pipeline-review` | Analyze pipeline health — prioritize deals, flag risks, get weekly action plan |

모든 명령은 note paste, CSV upload, 상황 설명만으로 **standalone** 동작하며, MCP connector를 연결하면 더 강력해집니다.

## 스킬

관련 상황에서 Claude가 자동으로 사용하는 domain knowledge입니다:

| 스킬 | 설명 |
|---|---|
| `account-research` | Company 또는 person을 조사하고 실행 가능한 sales intel을 얻습니다. Web search만으로 standalone 동작하며 enrichment tool 또는 CRM을 연결하면 더 강력해집니다. "research [company]", "look up [person]", "intel on [prospect]", "who is [name] at [company]", "tell me about [company]"에서 트리거됩니다. |
| `call-prep` | Account context, attendee research, suggested agenda로 sales call을 준비합니다. User input과 web research만으로 standalone 동작하며 CRM, email, chat, transcript를 연결하면 더 강력해집니다. "prep me for my call with [company]", "I'm meeting with [company] prep me", "call prep [company]"에서 트리거됩니다. |
| `daily-briefing` | Prioritized sales briefing으로 하루를 시작합니다. Meeting과 priority를 알려주면 standalone으로 동작하고 calendar, CRM, email을 연결하면 더 강력해집니다. "morning briefing", "daily brief", "what's on my plate today", "prep my day", "start my day"에서 트리거됩니다. |
| `draft-outreach` | Prospect를 조사한 뒤 개인화된 outreach 초안을 작성합니다. 기본적으로 web research만으로 동작하며 enrichment 및 CRM을 연결하면 더 강력해집니다. "draft outreach to [person/company]", "write cold email to [prospect]", "reach out to [name]"에서 트리거됩니다. |
| `competitive-intelligence` | Competitor를 조사하고 interactive battlecard를 만듭니다. Clickable competitor card와 comparison matrix가 있는 HTML artifact를 출력합니다. "competitive intel", "research competitors", "how do we compare to [competitor]", "battlecard for [competitor]", "what's new with [competitor]"에서 트리거됩니다. |
| `create-an-asset` | Deal context에서 tailored sales asset(landing page, deck, one-pager, workflow demo)을 생성합니다. Prospect, audience, goal을 설명하면 customer에게 공유할 ready-to-share polished branded asset을 얻습니다. |

## 예시 워크플로

### After a Call

```
/call-summary
```

Paste your notes or transcript. Get a structured summary, action items with owners, and a draft follow-up email. If CRM is connected, offers to log the activity and create tasks.

### Weekly Forecast

```
/forecast
```

Upload a CSV export from your CRM (or paste your deals). Tell me your quota and timeline. Get a weighted forecast with best/likely/worst scenarios, commit vs. upside breakdown, and gap analysis.

### Pipeline Review

```
/pipeline-review
```

Upload a CSV or describe your pipeline. Get a health score, deal prioritization, risk flags (stale deals, past close dates, single-threaded), and a weekly action plan.

### Researching a Prospect

Just ask naturally:
```
Research Acme Corp before my call tomorrow
```

The `account-research` skill triggers automatically and gives you company overview, key contacts, recent news, and recommended approach.

### Drafting Outreach

```
Draft an email to the VP of Engineering at TechStart
```

The `draft-outreach` skill researches the prospect first, then generates personalized outreach with multiple angles.

### Competitive Intel

```
How do we compare to Competitor X?
```

The `competitive-intelligence` skill researches both companies and builds a differentiation matrix with talk tracks.

## 단독 사용과 통합 사용

모든 명령과 스킬은 integration 없이도 동작합니다:

| 할 수 있는 일 | 단독 사용 | 통합 사용 시 |
|-----------------|------------|-------------------|
| Process call notes | Paste notes/transcript | Transcripts MCP (e.g. Gong, Fireflies) |
| Forecast pipeline | Upload CSV, paste deals | CRM MCP |
| Review pipeline | Upload CSV, describe deals | CRM MCP |
| Research prospects | Web search | Enrichment MCP (e.g. Clay, ZoomInfo) |
| Prep for calls | Describe meeting | CRM, Email, Calendar MCPs |
| Draft outreach | Web search + your context | CRM, Email MCPs |
| Competitive intel | Web search | CRM (win/loss data), Docs (battlecards) |

## MCP 통합

> If you see unfamiliar placeholders or need to check which tools are connected, see [CONNECTORS.md](CONNECTORS.md).

더 풍부한 경험을 위해 tool을 연결하세요:

| 범주 | 예시 | 가능해지는 일 |
|---|---|---|
| **CRM** | HubSpot, Close | Pipeline data, account history, contact records |
| **Transcripts** | Fireflies, Gong, Chorus | Call recordings, transcripts, key moments |
| **Enrichment** | Clay, ZoomInfo, Apollo | Company and contact data enrichment |
| **Chat** | Slack, Teams | Internal discussions, colleague intel |

Email, calendar, 추가 CRM option을 포함한 지원 integration 전체 목록은 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

## 설정

개인화를 위해 `settings.local.json` file을 만드세요:

- **Cowork**: Save it in any folder you've shared with Cowork (via the folder picker). The plugin finds it automatically.
- **Claude Code**: Save it at `sales/.claude/settings.local.json`.

```json
{
  "name": "Your Name",
  "title": "Account Executive",
  "company": "Your Company",
  "quota": {
    "annual": 1000000,
    "quarterly": 250000
  },
  "product": {
    "name": "Your Product",
    "value_props": [
      "Key value proposition 1",
      "Key value proposition 2"
    ],
    "competitors": [
      "Competitor A",
      "Competitor B"
    ]
  }
}
```

설정되어 있지 않으면 플러그인이 이 정보를 대화형으로 물어봅니다.
