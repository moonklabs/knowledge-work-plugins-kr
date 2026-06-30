# HR 플러그인

People operations plugin입니다. 주로 Anthropic의 agentic desktop application인 [Cowork](https://claude.com/product/cowork)용으로 설계되었지만 Claude Code에서도 동작합니다. Recruiting, onboarding, performance management, policy guidance, compensation analysis를 돕습니다. 어떤 HR team에서도 사용할 수 있으며, 입력만으로 standalone 동작하고 HRIS, ATS 및 다른 tool을 연결하면 더 강력해집니다.

## 설치

```bash
claude plugins add knowledge-work-plugins/human-resources
```

## 명령

Slash command로 호출하는 명시적 workflow입니다:

| 명령 | 설명 |
|---|---|
| `/draft-offer` | Draft an offer letter with compensation details, start date, and terms |
| `/onboarding` | Generate an onboarding checklist and first-week plan for a new hire |
| `/performance-review` | Structure a performance review — self-assessment prompts, manager template, calibration prep |
| `/policy-lookup` | Find and explain company policies — PTO, benefits, expense, travel, remote work |
| `/comp-analysis` | Analyze compensation data — benchmarking, band placement, equity refresh modeling |
| `/people-report` | Generate headcount, attrition, diversity, or org health reports |

All commands work **standalone** (provide context and details) and get **supercharged** with MCP connectors.

## 스킬

관련 상황에서 Claude가 자동으로 사용하는 domain knowledge입니다:

| 스킬 | 설명 |
|---|---|
| `recruiting-pipeline` | Recruiting pipeline stage를 track하고 manage합니다. "recruiting update", "candidate pipeline", "how many candidates", "hiring status" 또는 sourcing, screening, interviewing, offer extension 논의 시 트리거됩니다. |
| `employee-handbook` | Answer questions about company policies, benefits, and procedures |
| `compensation-benchmarking` | Benchmark compensation against market data — base, equity, total comp |
| `org-planning` | Headcount planning, org design, team structure optimization을 수행합니다. "org planning", "headcount plan", "team structure", "reorg", "who should we hire next" 또는 team size, reporting structure, organizational design을 고민할 때 트리거됩니다. |
| `people-analytics` | Analyze workforce data — attrition trends, engagement signals, diversity metrics |
| `interview-prep` | Competency-based question과 scorecard가 포함된 structured interview plan을 만듭니다. "interview plan for", "interview questions for", "how should we interview", "scorecard for" 또는 candidate interview 준비 요청에서 트리거됩니다. |

## 예시 워크플로

### Drafting an Offer

```
/draft-offer
```

Tell me the role, level, location, and comp details. Get a complete offer letter draft with terms, equity breakdown, and benefits summary.

### Onboarding a New Hire

```
/onboarding
```

Provide the new hire's name, role, team, and start date. Get a comprehensive onboarding checklist, first-week calendar, tool access list, and buddy assignment template.

### Preparing for Performance Reviews

```
/performance-review
```

Get templates for self-assessments, manager reviews, and calibration. I'll help structure feedback that's specific, actionable, and fair.

### Understanding Benefits

Just ask naturally:
```
What's our parental leave policy?
```

The `employee-handbook` skill triggers automatically and searches your connected knowledge base for the answer.

### Compensation Benchmarking

```
/comp-analysis
```

Upload comp data or describe your bands. Get market comparisons, band placement analysis, and recommendations for adjustments.

## 단독 사용과 통합 사용

모든 명령과 스킬은 integration 없이도 동작합니다:

| 할 수 있는 일 | 단독 사용 | 통합 사용 시 |
|-----------------|------------|-------------------|
| Draft offers | Provide details manually | HRIS, ATS for auto-fill |
| Onboarding checklists | Describe your process | HRIS, Knowledge base for templates |
| Performance reviews | Manual input | HRIS for review history |
| Policy lookup | Paste handbook content | Knowledge base |
| Comp analysis | Upload CSV, describe bands | Compensation data MCP |
| People reports | Upload data | HRIS for live data |

## MCP 통합

> If you see unfamiliar placeholders or need to check which tools are connected, see [CONNECTORS.md](CONNECTORS.md).

더 풍부한 경험을 위해 tool을 연결하세요:

| 범주 | 예시 | 가능해지는 일 |
|---|---|---|
| **HRIS** | Workday, BambooHR, Rippling | Employee data, org structure, PTO balances |
| **ATS** | Greenhouse, Lever, Ashby | Candidate pipeline, interview schedules, offer tracking |
| **Compensation** | Pave, Radford | Market benchmarks, comp band data |
| **Chat** | Slack, Teams | Team announcements, candidate coordination |
| **Calendar** | Google Calendar, Microsoft 365 | Interview scheduling, onboarding calendar |
| **Email** | Gmail, Microsoft 365 | Offer letters, candidate communications |

지원되는 integration 전체 목록은 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

## 설정

개인화를 위해 `settings.local.json` file을 만드세요:

- **Cowork**: Save it in any folder you've shared with Cowork (via the folder picker). The plugin finds it automatically.
- **Claude Code**: Save it at `human-resources/.claude/settings.local.json`.

```json
{
  "company": "Your Company",
  "headquarters": "City, State",
  "employeeCount": 500,
  "benefits": {
    "healthInsurance": "Provider Name",
    "pto": "Unlimited / X days",
    "parentalLeave": "X weeks"
  },
  "compensation": {
    "currency": "USD",
    "equityType": "RSU / Options",
    "vestingSchedule": "4 years, 1 year cliff"
  }
}
```

설정되어 있지 않으면 플러그인이 이 정보를 대화형으로 물어봅니다.
