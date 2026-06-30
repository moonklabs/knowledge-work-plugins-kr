# Sales 플러그인

영업 생산성 플러그인입니다. 주로 Anthropic의 agentic 데스크톱 애플리케이션인 [Cowork](https://claude.com/product/cowork)용으로 설계되었지만 Claude Code에서도 동작합니다. 잠재 고객 발굴, 아웃리치, 파이프라인 관리, 통화 준비, 거래 전략을 돕습니다. 어떤 영업 팀에서도 사용할 수 있으며, 웹 검색과 사용자 입력만으로 단독 동작하고 CRM, 이메일 및 다른 도구를 연결하면 더 강력해집니다.

## 설치

```bash
claude plugins add knowledge-work-plugins/sales
```

## 명령

슬래시 명령으로 호출하는 명시적 워크플로입니다.

| 명령 | 설명 |
|---|---|
| `/call-summary` | 통화 메모 또는 전사를 처리해 액션 아이템을 추출하고 후속 조치 초안과 내부 요약을 생성합니다 |
| `/forecast` | CSV를 업로드하거나 파이프라인을 설명하면 할당량 기준 가중 영업 전망을 생성합니다 |
| `/pipeline-review` | 파이프라인 상태를 분석해 거래 우선순위, 위험 표시, 주간 실행 계획을 제공합니다 |

모든 명령은 메모 붙여넣기, CSV 업로드, 상황 설명만으로 **단독** 동작하며, MCP connector를 연결하면 더 강력해집니다.

## 스킬

관련 상황에서 Claude가 자동으로 사용하는 도메인 지식입니다.

| 스킬 | 설명 |
|---|---|
| `account-research` | 회사 또는 사람을 조사하고 실행 가능한 영업 정보를 얻습니다. Web search만으로 단독 동작하며 enrichment tool 또는 CRM을 연결하면 더 강력해집니다. "research [company]", "look up [person]", "intel on [prospect]", "who is [name] at [company]", "tell me about [company]"에서 트리거됩니다. |
| `call-prep` | 계정 맥락, 참석자 조사, 추천 안건으로 영업 통화를 준비합니다. 사용자 입력과 웹 리서치만으로 단독 동작하며 CRM, email, chat, transcript를 연결하면 더 강력해집니다. "prep me for my call with [company]", "I'm meeting with [company] prep me", "call prep [company]"에서 트리거됩니다. |
| `daily-briefing` | 우선순위가 정해진 영업 브리핑으로 하루를 시작합니다. 미팅과 우선순위를 알려주면 단독으로 동작하고 calendar, CRM, email을 연결하면 더 강력해집니다. "morning briefing", "daily brief", "what's on my plate today", "prep my day", "start my day"에서 트리거됩니다. |
| `draft-outreach` | Prospect를 조사한 뒤 개인화된 outreach 초안을 작성합니다. 기본적으로 web research만으로 동작하며 enrichment 및 CRM을 연결하면 더 강력해집니다. "draft outreach to [person/company]", "write cold email to [prospect]", "reach out to [name]"에서 트리거됩니다. |
| `competitive-intelligence` | 경쟁사를 조사하고 인터랙티브 배틀카드를 만듭니다. 클릭 가능한 경쟁사 카드와 비교 매트릭스가 있는 HTML artifact를 출력합니다. "competitive intel", "research competitors", "how do we compare to [competitor]", "battlecard for [competitor]", "what's new with [competitor]"에서 트리거됩니다. |
| `create-an-asset` | 거래 맥락에서 맞춤형 영업 자료(landing page, deck, one-pager, workflow demo)를 생성합니다. Prospect, audience, goal을 설명하면 고객에게 공유할 수 있는 다듬어진 브랜드 자료를 얻습니다. |

## 예시 워크플로

### 통화 이후

```
/call-summary
```

메모 또는 transcript를 붙여넣으면 구조화된 요약, 담당자가 있는 액션 아이템, 후속 email 초안을 받습니다. CRM이 연결되어 있으면 activity logging과 task 생성을 제안합니다.

### 주간 전망

```
/forecast
```

CRM에서 내보낸 CSV를 업로드하거나 deal을 붙여넣습니다. 할당량과 일정을 알려주면 best/likely/worst scenario, commit vs. upside breakdown, gap analysis가 포함된 weighted forecast를 받습니다.

### 파이프라인 검토

```
/pipeline-review
```

CSV를 업로드하거나 pipeline을 설명합니다. Health score, deal prioritization, stale deal/past close date/single-threaded 같은 위험 표시와 주간 실행 계획을 받습니다.

### Prospect 조사

자연스럽게 요청하세요.
```
내일 call 전에 Acme Corp 조사해줘
```

`account-research` skill이 자동으로 trigger되어 company overview, key contact, recent news, recommended approach를 제공합니다.

### 아웃리치 작성

```
TechStart VP of Engineering에게 보낼 email 초안 작성해줘
```

`draft-outreach` skill이 prospect를 먼저 조사한 뒤 여러 angle의 personalized outreach를 생성합니다.

### 경쟁 정보

```
Competitor X와 비교하면 우리는 어때?
```

`competitive-intelligence` skill이 두 회사를 조사하고 talk track이 포함된 differentiation matrix를 만듭니다.

## 단독 사용과 통합 사용

모든 명령과 스킬은 통합 도구 없이도 동작합니다.

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

> 익숙하지 않은 placeholder가 보이거나 연결된 도구를 확인해야 한다면 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

더 풍부한 경험을 위해 도구를 연결하세요.

| 범주 | 예시 | 가능해지는 일 |
|---|---|---|
| **CRM** | HubSpot, Close | Pipeline data, account history, contact records |
| **Transcripts** | Fireflies, Gong, Chorus | Call recordings, transcripts, key moments |
| **Enrichment** | Clay, ZoomInfo, Apollo | Company and contact data enrichment |
| **Chat** | Slack, Teams | Internal discussions, colleague intel |

Email, calendar, 추가 CRM option을 포함한 지원 통합 전체 목록은 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

## 설정

개인화를 위해 `settings.local.json` 설정 파일을 만드세요.

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
