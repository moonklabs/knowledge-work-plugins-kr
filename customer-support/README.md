# Customer Support 플러그인

Support team을 위한 customer support plugin입니다. 주로 Anthropic의 agentic desktop application인 [Cowork](https://claude.com/product/cowork)용으로 설계되었지만 Claude Code에서도 동작합니다. Ticket triage, escalation management, response drafting, customer research, knowledge base authoring을 제공합니다.

## 설치

```
claude plugins add knowledge-work-plugins/customer-support
```

## 주요 기능

이 plugin은 Claude를 customer support co-pilot으로 바꿉니다. 다음을 도와줍니다:

- **Incoming ticket triage**: structured categorization, priority assessment, routing recommendation을 제공합니다.
- **Customer question research**: 여러 source의 정보를 confidence scoring과 함께 종합합니다.
- **Professional response draft**: 상황, 긴급도, communication channel에 맞춘 답변을 작성합니다.
- **Escalation packaging**: engineering/product team을 위해 full context, reproduction step, business impact를 정리합니다.
- **KB article 작성**: 해결된 issue를 문서화해 future ticket volume을 줄입니다.

## 명령

| 명령 | 설명 |
|---|---|
| `/triage` | Support ticket 또는 customer issue를 categorize, prioritize, route합니다. |
| `/research` | Customer question 또는 topic에 대해 multi-source research를 수행합니다. |
| `/draft-response` | 어떤 상황에서도 customer-facing response 초안을 작성합니다. |
| `/escalate` | Engineering, product, leadership용 escalation을 package합니다. |
| `/kb-article` | Resolved issue에서 knowledge base article 초안을 작성합니다. |

## 스킬

| 스킬 | 설명 |
|---|---|
| `ticket-triage` | Support ticket 또는 customer issue를 triage하고 prioritize합니다. New ticket이 들어와 categorization, P1-P4 priority assignment, 담당 team 결정, routing 전 duplicate/known issue 확인이 필요할 때 사용합니다. |
| `customer-research` | Source attribution이 포함된 customer question/topic에 대한 multi-source research입니다. Customer가 lookup이 필요한 질문을 하거나, bug가 이전에 report되었는지 조사하거나, specific account에 이전에 무엇을 안내했는지 확인하거나, response 초안 전 background를 수집할 때 사용합니다. |
| `response-drafting` | Communication best practice, tone guideline, common scenario template을 제공합니다. |
| `escalation` | Escalation tier, structured escalation format, impact assessment, follow-up cadence를 다룹니다. |
| `knowledge-management` | Article structure standard, searchability를 위한 writing, review cadence, maintenance를 다룹니다. |

## Data Source

> If you see unfamiliar placeholders or need to check which tools are connected, see [CONNECTORS.md](CONNECTORS.md).

Connect your support platform, knowledge base, and communication tools for the best experience. Without them, provide customer context manually.

**Included MCP connections:**
- Chat (Slack) for internal discussions and customer channel context
- Support platform (Intercom) for ticket history and customer conversations
- CRM (HubSpot) for account details and contact information
- Knowledge base (Guru, Notion) for internal documentation and runbooks
- Project tracker (Atlassian) for bug reports and feature requests
- Email and cloud storage (Microsoft 365) for customer correspondence

**Additional options:**
- See [CONNECTORS.md](CONNECTORS.md) for alternative tools in each category

## 예시 워크플로

### Triaging an Incoming Ticket

```
You: /triage Customer says their dashboard has been showing a blank
     page since this morning — they're on the Enterprise plan

Claude: Analyzing the issue...

## Triage: Dashboard blank page — Enterprise customer

**Category:** Bug
**Priority:** P2 — High (core feature broken, enterprise customer)
**Product area:** Dashboard / Frontend

### Routing Recommendation
Route to: Tier 2 Support (investigate browser/account specifics,
then escalate to Engineering if reproducible)

### Suggested Initial Response
[Draft acknowledging the issue, asking clarifying questions
about browser and any recent changes, offering status page link]
```

### Researching a Customer Question

```
You: /research Does our platform support SSO with Okta?

Claude: Let me check across available sources...

[Searches docs, knowledge base, previous conversations,
 and web sources — synthesizes answer with confidence
 level and source attribution]
```

### Packaging an Escalation

```
You: /escalate API returning 500 errors intermittently — 3
     Enterprise customers reported this week

Claude: I'll gather context and build the escalation brief...

[Pulls ticket history, identifies affected accounts,
 structures reproduction steps, assesses business impact,
 and generates a structured escalation for engineering]
```

### Drafting a Response

```
You: /draft-response Customer escalation — their integration
     has been down for 2 days and they're frustrated

Claude: [Generates empathetic, professional response acknowledging
 the frustration, taking ownership, providing a clear action
 plan with timeline, and offering a call]
```

### Writing a KB Article

```
You: /kb-article How to configure webhook notifications —
     just resolved this for the third customer this month

Claude: [Generates a structured how-to article with prerequisites,
 step-by-step instructions, verification steps, and common
 issues — optimized for search]
```

## 설정

The plugin works out of the box with the included MCP connections. For the richest experience, connect additional data sources through your Claude settings:

1. **Support platform**: Add your ticketing system for ticket history and customer context
2. **Knowledge base**: Add your wiki for internal documentation and existing KB articles
3. **Project tracker**: Add your issue tracker for bug reports and feature requests
4. **CRM**: Add your CRM for account details and contact information

Without these connections, the plugin will ask you to provide context manually and offer frameworks and templates you can fill in with your own data.
