# Customer Support 플러그인

Support team을 위한 customer support 플러그인입니다. 주로 Anthropic의 agentic 데스크톱 애플리케이션인 [Cowork](https://claude.com/product/cowork)용으로 설계되었지만 Claude Code에서도 동작합니다. 티켓 분류, 에스컬레이션 관리, 답변 초안 작성, 고객 조사, 지식 베이스 문서 작성을 제공합니다.

## 설치

```text
claude plugins add knowledge-work-plugins/customer-support
```

## 주요 기능

이 플러그인은 Claude를 customer support co-pilot으로 바꿉니다. 다음을 도와줍니다.

- **수신 티켓 분류**: 구조화된 분류, 우선순위 평가, 라우팅 추천을 제공합니다.
- **고객 질문 조사**: 여러 출처의 정보를 confidence scoring과 함께 종합합니다.
- **전문적인 답변 초안**: 상황, 긴급도, 커뮤니케이션 채널에 맞춘 답변을 작성합니다.
- **Escalation package**: engineering/product team을 위해 전체 맥락, 재현 단계, 비즈니스 영향을 정리합니다.
- **KB article 작성**: 해결된 이슈를 문서화해 향후 티켓 양을 줄입니다.

## 명령

| 명령 | 설명 |
|---|---|
| `/triage` | Support ticket 또는 customer issue를 분류하고 우선순위화한 뒤 라우팅합니다. |
| `/research` | Customer question 또는 topic에 대해 여러 출처를 조사합니다. |
| `/draft-response` | 어떤 상황에서도 customer-facing response 초안을 작성합니다. |
| `/escalate` | Engineering, product, leadership용 escalation을 package합니다. |
| `/kb-article` | 해결된 이슈에서 knowledge base article 초안을 작성합니다. |

## 스킬

| 스킬 | 설명 |
|---|---|
| `ticket-triage` | Support ticket 또는 customer issue를 triage하고 prioritize합니다. 새 티켓이 들어와 분류, P1-P4 우선순위 지정, 담당 team 결정, routing 전 duplicate/known issue 확인이 필요할 때 사용합니다. |
| `customer-research` | 출처 표시가 포함된 customer question/topic에 대한 여러 출처 조사입니다. 고객 질문에 조사가 필요하거나, bug가 이전에 report되었는지 확인하거나, 특정 account에 이전에 무엇을 안내했는지 확인하거나, response 초안 전 배경을 수집할 때 사용합니다. |
| `response-drafting` | 커뮤니케이션 best practice, tone guideline, common scenario template을 제공합니다. |
| `escalation` | Escalation tier, 구조화된 escalation 형식, 영향 평가, follow-up cadence를 다룹니다. |
| `knowledge-management` | Article structure standard, 검색 가능성을 위한 writing, review cadence, maintenance를 다룹니다. |

## 데이터 출처

> 익숙하지 않은 placeholder가 보이거나 연결된 도구를 확인해야 한다면 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

가장 좋은 경험을 위해 support platform, knowledge base, 커뮤니케이션 도구를 연결하세요. 연결하지 않았다면 고객 맥락을 수동으로 제공하면 됩니다.

**포함된 MCP 연결:**

- 채팅(Slack): 내부 논의 및 customer channel context
- Support platform(Intercom): ticket history 및 customer conversation
- CRM(HubSpot): account detail 및 contact information
- Knowledge base(Guru, Notion): internal documentation 및 runbook
- Project tracker(Atlassian): bug report 및 feature request
- Email 및 cloud storage(Microsoft 365): customer correspondence

**추가 옵션:**

- 범주별 대체 도구는 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

## 예시 워크플로

### 수신 티켓 분류

```text
You: /triage Enterprise 플랜 고객인데 오늘 아침부터 dashboard가 빈 페이지로 보인다고 합니다

Claude: 이슈를 분석합니다...

## Triage: Dashboard blank page — Enterprise customer

**Category:** Bug
**Priority:** P2 — High (core feature broken, enterprise customer)
**Product area:** Dashboard / Frontend

### Routing Recommendation
Route to: Tier 2 Support (브라우저/계정 특이사항을 조사하고,
재현 가능하면 Engineering으로 escalate)

### Suggested Initial Response
[이슈를 인지했음을 알리고, 브라우저와 최근 변경 사항에 대한
확인 질문을 하며, status page link를 제공하는 초안]
```

### Customer question 조사

```text
You: /research 우리 platform은 Okta SSO를 지원하나요?

Claude: 사용 가능한 출처를 확인해볼게요...

[문서, knowledge base, 이전 대화, 웹 출처를 검색하고
 confidence level과 출처 표시가 포함된 답변으로 종합]
```

### Escalation package 작성

```text
You: /escalate API가 간헐적으로 500 error를 반환합니다 — 이번 주 Enterprise 고객 3곳이 보고했습니다

Claude: 맥락을 모아 escalation brief를 만들겠습니다...

[ticket history를 가져오고, 영향을 받은 account를 식별하고,
 재현 단계를 구조화하며, business impact를 평가한 뒤
 engineering용 structured escalation을 생성]
```

### Response 초안 작성

```text
You: /draft-response 고객 escalation — integration이 2일째 down되어 고객이 frustrated 상태입니다

Claude: [불편을 인정하고 ownership을 보이며,
 timeline이 포함된 명확한 action plan과 call 제안을 담은
 공감 어린 전문 답변을 생성]
```

### KB article 작성

```text
You: /kb-article Webhook notification 설정 방법 — 이번 달 세 번째 고객에게 방금 해결해줬습니다

Claude: [전제 조건, 단계별 안내,
 verification step, common issue가 포함된 구조화된 how-to article을 생성하고
 search에 최적화]
```

## 설정

플러그인은 포함된 MCP 연결만으로 바로 동작합니다. 더 풍부한 경험을 위해 Claude settings에서 추가 데이터 출처를 연결하세요.

1. **Support platform**: Ticket history와 고객 맥락을 위해 ticketing system을 추가합니다.
2. **Knowledge base**: Internal documentation과 existing KB article을 위해 wiki를 추가합니다.
3. **Project tracker**: Bug report와 feature request를 위해 issue tracker를 추가합니다.
4. **CRM**: Account detail과 contact information을 위해 CRM을 추가합니다.

이 연결이 없으면 플러그인은 context를 수동으로 제공하도록 요청하고, 사용자가 own data로 채울 수 있는 framework와 template을 제공합니다.
