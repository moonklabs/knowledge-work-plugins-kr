# Sales 플러그인

Sales productivity 플러그인입니다. 주로 Anthropic의 agentic 데스크톱 애플리케이션인 [Cowork](https://claude.com/product/cowork)용으로 설계되었지만 Claude Code에서도 동작합니다. Prospecting, outreach, pipeline management, call preparation, deal strategy를 돕습니다. 어떤 sales team에서도 사용할 수 있으며, web search와 사용자 입력만으로 단독 동작하고 CRM, email 및 다른 tool을 연결하면 더 강력해집니다.

## 설치

```bash
claude plugins add knowledge-work-plugins/sales
```

## 명령

Slash command로 호출하는 명시적 워크플로입니다:

| 명령 | 설명 |
|---|---|
| `/call-summary` | Call note 또는 transcript를 처리해 action item을 추출하고 follow-up 초안과 internal summary를 생성합니다 |
| `/forecast` | CSV를 upload하거나 pipeline을 설명하면 quota 기준 weighted sales forecast를 생성합니다 |
| `/pipeline-review` | Pipeline health를 분석해 deal 우선순위, risk flag, weekly action plan을 제공합니다 |

모든 명령은 note paste, CSV upload, 상황 설명만으로 **단독** 동작하며, MCP connector를 연결하면 더 강력해집니다.

## 스킬

관련 상황에서 Claude가 자동으로 사용하는 domain knowledge입니다:

| 스킬 | 설명 |
|---|---|
| `account-research` | Company 또는 person을 조사하고 실행 가능한 sales intel을 얻습니다. Web search만으로 standalone 동작하며 enrichment tool 또는 CRM을 연결하면 더 강력해집니다. "research [company]", "look up [person]", "intel on [prospect]", "who is [name] at [company]", "tell me about [company]"에서 트리거됩니다. |
| `call-prep` | Account context, attendee research, suggested agenda로 sales call을 준비합니다. User input과 web research만으로 단독 동작하며 CRM, email, chat, transcript를 연결하면 더 강력해집니다. "prep me for my call with [company]", "I'm meeting with [company] prep me", "call prep [company]"에서 트리거됩니다. |
| `daily-briefing` | Prioritized sales briefing으로 하루를 시작합니다. Meeting과 priority를 알려주면 standalone으로 동작하고 calendar, CRM, email을 연결하면 더 강력해집니다. "morning briefing", "daily brief", "what's on my plate today", "prep my day", "start my day"에서 트리거됩니다. |
| `draft-outreach` | Prospect를 조사한 뒤 개인화된 outreach 초안을 작성합니다. 기본적으로 web research만으로 동작하며 enrichment 및 CRM을 연결하면 더 강력해집니다. "draft outreach to [person/company]", "write cold email to [prospect]", "reach out to [name]"에서 트리거됩니다. |
| `competitive-intelligence` | Competitor를 조사하고 interactive battlecard를 만듭니다. Clickable competitor card와 comparison matrix가 있는 HTML artifact를 출력합니다. "competitive intel", "research competitors", "how do we compare to [competitor]", "battlecard for [competitor]", "what's new with [competitor]"에서 트리거됩니다. |
| `create-an-asset` | Deal context에서 tailored sales asset(landing page, deck, one-pager, workflow demo)을 생성합니다. Prospect, audience, goal을 설명하면 customer에게 공유할 ready-to-share polished branded asset을 얻습니다. |

## 예시 워크플로

### Call 이후

```
/call-summary
```

Note 또는 transcript를 붙여넣으면 구조화된 summary, owner가 있는 action item, follow-up email draft를 받습니다. CRM이 연결되어 있으면 activity logging과 task 생성을 제안합니다.

### Weekly forecast

```
/forecast
```

CRM에서 export한 CSV를 upload하거나 deal을 paste합니다. Quota와 timeline을 알려주면 best/likely/worst scenario, commit vs. upside breakdown, gap analysis가 포함된 weighted forecast를 받습니다.

### Pipeline review

```
/pipeline-review
```

CSV를 upload하거나 pipeline을 설명합니다. Health score, deal prioritization, stale deal/past close date/single-threaded 같은 risk flag, weekly action plan을 받습니다.

### Prospect 조사

자연스럽게 요청하세요.
```
내일 call 전에 Acme Corp 조사해줘
```

`account-research` skill이 자동으로 trigger되어 company overview, key contact, recent news, recommended approach를 제공합니다.

### Outreach 작성

```
TechStart VP of Engineering에게 보낼 email 초안 작성해줘
```

`draft-outreach` skill이 prospect를 먼저 조사한 뒤 여러 angle의 personalized outreach를 생성합니다.

### Competitive intel

```
Competitor X와 비교하면 우리는 어때?
```

`competitive-intelligence` skill이 두 회사를 조사하고 talk track이 포함된 differentiation matrix를 만듭니다.

## 단독 사용과 통합 사용

모든 명령과 스킬은 integration 없이도 동작합니다:

| 할 수 있는 일 | 단독 사용 | 통합 사용 시 |
|-----------------|------------|-------------------|
| Call note 처리 | Note/transcript 붙여넣기 | Transcripts MCP(예: Gong, Fireflies) |
| Pipeline forecast | CSV upload, deal paste | CRM MCP |
| Pipeline review | CSV upload, deal 설명 | CRM MCP |
| Prospect 조사 | Web search | Enrichment MCP(예: Clay, ZoomInfo) |
| Call prep | Meeting 설명 | CRM, Email, Calendar MCP |
| Outreach 초안 | Web search + 사용자 context | CRM, Email MCP |
| Competitive intel | Web search | CRM(win/loss data), Docs(battlecards) |

## MCP 통합

> 익숙하지 않은 placeholder가 보이거나 연결된 tool을 확인해야 한다면 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

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

- **Cowork**: Folder picker로 Cowork와 공유한 아무 folder에 저장합니다. Plugin이 자동으로 찾습니다.
- **Claude Code**: `sales/.claude/settings.local.json`에 저장합니다.

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
