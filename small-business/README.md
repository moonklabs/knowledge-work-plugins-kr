# Small Business 플러그인

[Cowork](https://claude.com/product/cowork)용으로 미리 구성된 small business workflow 모음입니다. Cowork는 Anthropic의 agentic desktop application이며, 이 플러그인은 Claude Code에서도 동작합니다. 한 번 설치하면 15개의 building-block skill, 15개의 바로 쓸 수 있는 workflow, 자연어 요청을 이해하는 router를 사용할 수 있습니다.

무엇을 외울 필요가 없습니다. Claude에게 필요한 일을 그대로 말하면 됩니다. 예: "급여를 지급할 수 있을지 걱정돼", "고객이 화가 났어", "가격을 얼마로 잡아야 할까?". Claude가 적절한 workflow를 찾고 단계별로 안내합니다. 모든 workflow는 action 전에 멈추므로, owner의 승인 없이 실행되는 일은 없습니다.

> **중요**: 이 플러그인은 small business workflow를 보조하지만 financial, tax, legal, HR advice를 제공하지 않습니다. 모든 output은 사용 전 직접 검토해야 하며, 필요한 경우 자격 있는 전문가의 검토를 받아야 합니다.

## 설치

### Cowork

[claude.com/plugins](https://claude.com/plugins/)에서 설치합니다.

### Claude Code

```bash
claude plugin marketplace add anthropics/knowledge-work-plugins
claude plugin install small-business@knowledge-work-plugins
```

설치가 끝나면 **"set me up"**이라고 말해 `smb-onboard` skill을 실행하세요. Claude가 사업, pain point, 이미 사용하는 tool을 이해하도록 도와줍니다.

## 연결하면 좋은 도구

`/smb-onboard`를 실행하거나 Claude에게 "set me up"이라고 말하세요.

**핵심 도구** (가장 좋은 경험을 위해 먼저 연결):
- **QuickBooks** — cash forecasts, margins, month-end close, tax prep 등 financial workflow를 구동합니다
- **PayPal** — transaction data, invoices, disputes, refunds
- **HubSpot** — CRM, leads, campaigns, customer support tickets

**Marketing & communication:**
- **Canva** — on-brand social 및 email asset 생성
- **Gmail / Outlook** — email draft, ticket handling, contract review
- **Google Calendar / Outlook Calendar** — meeting prep, call blocking, weekly commitments
- **Slack** — brief delivery 및 notifications

**선택 도구** (연결하면 더 깊은 분석 가능):
- **Stripe** — payment and subscription data
- **Square** — POS transaction data
- **Google Drive / OneDrive** — file storage and templates
- **DocuSign** — contract review from pending envelopes

처음부터 전부 연결할 필요는 없습니다. 한두 개만 연결해도 바로 가치가 보이며, 다른 tool을 연결하면 더 많은 기능이 열릴 때 plugin이 알려줍니다.

## 작동 방식

세 계층이 함께 동작합니다.

1. **Skills** — building block입니다. 각 skill은 cash forecast, lead scoring, invoice reminder draft처럼 한 가지 일을 잘 처리합니다. 총 15개가 있습니다.

2. **Commands** — workflow입니다. Command는 여러 skill을 multi-step recipe로 연결하고, 실행 전 승인 checkpoint를 둡니다. 총 15개가 있습니다.

3. **Router** — front door입니다. Claude에게 자연어로 말하면 router가 적합한 workflow를 찾아 안내합니다. Command name을 외울 필요가 없습니다.

## 전체 15개 명령

Command는 여러 skill을 연결하는 workflow입니다. 각 command는 action 전에 approval checkpoint에서 멈춥니다.

### 자금 및 재무

| 명령 | 수행 내용 | 이렇게 말해보세요 | 사용 skill | 필수 | 선택 |
|---|---|---|---|---|---|
| `/plan-payroll` | Cash forecast + overdue invoice chase so you know payroll is covered. | "can I make payroll", "cash is tight", "who owes me money" | cash-flow-snapshot, invoice-chase | QuickBooks | PayPal, Stripe, Square, Mail |
| `/month-heads-up` | 30-day cash outlook with early risk flags. | "what does next month look like", "cash forecast", "runway" | cash-flow-snapshot | QuickBooks | PayPal |
| `/close-month` | Month-end close: reconcile, flag gaps, write P&L, export packet. | "close the books", "month-end", "reconcile" | month-end-prep | QuickBooks | PayPal, Stripe, Square |
| `/price-check` | Margin-by-product table and three pricing scenarios. | "what are my margins", "should I raise prices", "cost per unit" | margin-analyzer | QuickBooks | PayPal |
| `/tax-prep` | Tax prep materials for your accountant (quarterly estimates or year-end 1099s). | "tax stuff", "estimated taxes", "1099s", "accountant needs..." | tax-season-organizer | QuickBooks | PayPal, Stripe |

### 영업 및 마케팅

| 명령 | 수행 내용 | 이렇게 말해보세요 | 사용 skill | 필수 | 선택 |
|---|---|---|---|---|---|
| `/call-list` | Top 5 leads to call today with talking points and calendar blocks. | "who should I call", "any hot leads", "pipeline" | lead-triage | HubSpot | Mail, Google Calendar |
| `/run-campaign` | End-to-end campaign: sales analysis → content brief → Canva assets → HubSpot send. | "run a campaign", "sales are down", "I need more customers" | content-strategy, canva-creator, lead-triage | HubSpot, Canva | QuickBooks, PayPal |
| `/sales-brief` | Top and bottom sellers with a 2-week content brief. | "what's selling", "what should I promote" | content-strategy | QuickBooks or PayPal | HubSpot |

### 고객 및 운영

| 명령 | 수행 내용 | 이렇게 말해보세요 | 사용 skill | 필수 | 선택 |
|---|---|---|---|---|---|
| `/customer-pulse-check` | Customer feedback themes with response templates. | "what are customers saying", "complaints", "reviews" | customer-pulse, ticket-deflector | PayPal or HubSpot | -- |
| `/handle-complaint` | End-to-end complaint resolution: pull context, draft response, suggest operational fix. | "a customer is upset", "handle this complaint", "angry email" | ticket-deflector, customer-pulse | -- (works with pasted text) | Gmail, HubSpot, PayPal |
| `/crm-cleanup` | HubSpot hygiene: stale deals, duplicates, missing fields — fixes what you approve. | "clean up the CRM", "HubSpot is a mess", "stale deals" | crm-maintenance | HubSpot | -- |
| `/review-contract` | Plain-English contract review with red flags and severity ratings. | "review this contract", "NDA", "should I sign this" | contract-review | -- (works with file upload) | DocuSign |

### Business intelligence

| 명령 | 수행 내용 | 이렇게 말해보세요 | 사용 skill | 필수 | 선택 |
|---|---|---|---|---|---|
| `/monday-brief` | Monday morning briefing: cash, sales, pipeline, week ahead, top 3 to-dos. | "Monday brief", "what's on my plate", "start of week" | business-pulse | -- (degrades gracefully) | QuickBooks, PayPal, HubSpot, Calendar, Gmail, Slack |
| `/friday-brief` | Friday end-of-week pulse: revenue vs last week, wins, and things to watch. | "end of week", "how'd we do", "Friday recap" | business-pulse | PayPal or HubSpot | -- |
| `/quarterly-review` | Full QBR narrative: revenue, margin, customer health, opportunities, risks. | "quarterly review", "board deck", "QBR" | business-pulse | QuickBooks | PayPal, HubSpot |

## 전체 15개 스킬

Skill은 atomic building block입니다. 각 skill은 한 가지 일을 잘 처리합니다.

### 자금 및 재무

| 스킬 | 수행 내용 | 이렇게 말해보세요 | 필수 | 선택 |
|---|---|---|---|---|
| **cash-flow-snapshot** | 30/60/90-day cash forecast with confidence bands and named risk flags. Chat summary + XLSX. | "forecast my cash flow", "will I make payroll", "runway", "cash crunch" | QuickBooks, PayPal, Stripe, or Square (any one) | Others as secondary sources |
| **invoice-chase** | Drafts overdue-invoice reminders matched to each customer's payment history and tone. Sends via PayPal with approval. | "who owes me money", "overdue invoices", "follow up on unpaid" | QuickBooks | PayPal, Stripe, Gmail |
| **margin-analyzer** | Unit economics by product or service with inflation benchmarks and three pricing scenarios. | "what are my margins", "should I raise prices", "costs eating into profit", "what to charge" | QuickBooks | PayPal, Square, CSV upload |
| **month-end-prep** | Month-end close: reconciles QB against payment processors, flags gaps, writes P&L narrative, exports close packet. | "close the month", "reconcile", "P&L", "why revenue changed" | QuickBooks | PayPal, Stripe, Square |
| **tax-season-organizer** | Quarterly estimated tax calc or year-end 1099-NEC prep with accountant handoff packet. | "quarterly taxes", "estimated tax payment", "1099s", "1099-NEC", "year-end tax prep" | QuickBooks | PayPal, Stripe |

### 영업 및 마케팅

| 스킬 | 수행 내용 | 이렇게 말해보세요 | 필수 | 선택 |
|---|---|---|---|---|
| **lead-triage** | Scores HubSpot leads by engagement, fit, and urgency to produce a ranked call list with talking points. | "prioritize leads", "who to call first", "pipeline" | HubSpot | Gmail, Google Calendar |
| **content-strategy** | Analyzes sales data to find top performers and slow movers, produces a prioritized 30-day content brief. | "what should I post", "content plan", "what's selling", "what to promote" | QuickBooks or PayPal | Square |
| **canva-creator** | Takes a content brief and executes the full campaign: posting calendar, Canva assets, caption copy, HubSpot staging. | "make the content", "generate the posts", "create the assets", "turn this into a campaign" | Canva, HubSpot | -- |

### 고객 및 운영

| 스킬 | 수행 내용 | 이렇게 말해보세요 | 필수 | 선택 |
|---|---|---|---|---|
| **customer-pulse** | Aggregates disputes, tickets, email sentiment, and reviews into a themes report with a "do these three things this week" list. | "how are customers feeling", "what people are saying", "disputes", "review analysis" | -- (degrades gracefully) | PayPal, HubSpot, Gmail |
| **ticket-deflector** | Reads a customer email or ticket, pulls order/refund status, drafts a tone-matched reply. Can issue PayPal refunds with approval. | "draft a response", "answer this customer", "where's my order", "I want a refund" | PayPal, HubSpot, Mail | Intercom, Square |
| **crm-maintenance** | Keeps HubSpot current: creates/updates contacts and deals, logs calls and notes, flags stale records. | "update the CRM", "log a call", "clean up HubSpot", "add context to a deal" | HubSpot | Gmail, Google Calendar |
| **contract-review** | Plain-English contract review with risk flags, severity ratings, and a marked-up redline DOCX. | "review this contract", "what am I signing", "flag any concerns", "check the payment terms" | -- (works with file upload) | Gmail, DocuSign |

### 채용

| 스킬 | 수행 내용 | 이렇게 말해보세요 | 필수 | 선택 |
|---|---|---|---|---|
| **job-post-builder** | Builds a complete hiring packet: job post, structured interview guide with scoring rubric, and offer letter template. | "help me hire", "write a job post", "job description", "open role", "interview questions", "draft an offer letter" | -- (works standalone) | DocuSign, Google Drive |

### Business intelligence 및 onboarding

| 스킬 | 수행 내용 | 이렇게 말해보세요 | 필수 | 선택 |
|---|---|---|---|---|
| **business-pulse** | One-page business snapshot: cash, sales, pipeline, commitments, watch-list, and the single most important thing needing attention today. | "how's the business doing", "snapshot", "weekly summary", "catch me up" | -- (degrades gracefully) | QuickBooks, PayPal, HubSpot, Google Calendar, Gmail, Slack |
| **smb-onboard** | Walks you through connecting tools, runs a demo recipe, captures your business context, and sets a weekly check-in cadence. | "set me up", "setup", "get started", "help me get set up", "I'm new to this", "what can you do" | -- | All connectors |

## 커스터마이징

이 workflow들은 generic starting point입니다. 실제 사업 방식에 맞게 customize하면 훨씬 유용해집니다.

- **Business context 추가** — Industry, products, customers, process를 skill file에 넣어 Claude가 맥락을 이해하게 합니다.
- **Threshold 조정** — `business-pulse`와 `cash-flow-snapshot`의 alert threshold를 규모에 맞게 조정합니다.
- **Connector 교체** — 실제 사용하는 tool을 skill이 바라보게 설정합니다.
