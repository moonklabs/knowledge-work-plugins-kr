# Common Room 플러그인

Common Room 기반 GTM 워크플로입니다. Account research, contact research, call prep, personalized outreach, prospecting, weekly briefing을 지원합니다.

## 개요

이 plugin은 Claude를 Common Room MCP server에 연결하고, 가장 일반적인 rep workflow를 다루는 여섯 가지 skill을 제공합니다. 모든 output은 실제 Common Room signal data에 근거합니다. 1st-party product signal, 2nd-party community signal, 3rd-party intent signal, RoomieAI와 Spark의 enrichment를 활용합니다.

## 요구 사항

- **Common Room MCP**(`mcp.commonroom.io/mcp`)가 연결되고 인증되어 있어야 합니다. 모든 plugin 기능의 기본 데이터 소스입니다.
- **Calendar connector**(선택) — `call-prep`와 `weekly-prep-brief`에서 자동 meeting lookup을 가능하게 합니다. 연결하지 않으면 두 skill 모두 meeting detail을 사용자에게 묻습니다.

## 스킬

Skill은 대화로 trigger됩니다. 원하는 일을 설명하면 Claude가 적절한 skill을 자동으로 불러옵니다.

| 스킬 | Trigger phrase |
|-------|----------------|
| `account-research` | Common Room 데이터를 사용해 회사를 조사합니다. 'research [company]', 'tell me about [domain]', 'pull up signals for [account]', 'what's going on with [company]' 또는 계정 수준 질문에서 트리거됩니다. |
| `contact-research` | Common Room 데이터를 사용해 특정 사람을 조사합니다. 'who is [name]', 'look up [email]', 'research [contact]', 'is [name] a warm lead' 또는 contact-level question에서 트리거됩니다. |
| `call-prep` | Common Room signal을 사용해 customer 또는 prospect call을 준비합니다. 'prep me for my call with [company]', 'prepare for a meeting with [company]', 'what should I know before talking to [company]' 또는 call preparation request에서 트리거됩니다. |
| `compose-outreach` | Common Room signal을 사용해 개인화된 outreach message를 생성합니다. 'draft outreach to [person]', 'write an email to [name]', 'compose a message for [contact]' 또는 outreach drafting request에서 트리거됩니다. |
| `prospect` | Common Room Prospector를 사용해 targeted account 또는 contact list를 만듭니다. 'find companies that match [criteria]', 'build a prospect list', 'find contacts at [type of company]', 'show me companies hiring [role]' 또는 list-building request에서 트리거됩니다. |
| `weekly-prep-brief` | 다음 7일 동안의 모든 external call에 대한 종합 weekly briefing을 생성합니다. 'weekly prep brief', 'prepare my week', 'what calls do I have this week', 'Monday prep' 또는 weekly planning request에서 트리거됩니다. |

## 명령

명시적 호출이 더 적합한 complex workflow용 command 두 가지입니다.

| 명령 | 사용법 |
|---------|-------|
| `/generate-account-plan <company>` | Stakeholder mapping, engagement analysis, opportunity, risk, action item이 포함된 종합 strategic account plan |
| `/weekly-brief [date range]` | Full weekly prep briefing을 생성합니다. 기본값은 다음 7일입니다 |

## 각 스킬의 output

**Account Research** — Full overview, targeted field question, honest sparse-data response, MCP data + LLM reasoning 결합 등 네 pattern을 처리합니다. Recent news용 web search를 포함하며 "My Segments"로 자동 scope합니다.

**Contact Research** — Email, name+company, social handle로 lookup합니다. Enriched identity, CRM field, score, website visit, activity history, Spark analysis, conversation starter를 반환합니다.

**Call Prep** — Company snapshot, attendee별 profile, signal highlight, tailored talking point, likely objection, recommended call outcome을 제공합니다. Gong/call recording activity를 우선하며, 연결된 경우 calendar-aware하게 동작합니다.

**Compose Outreach** — Common Room signal과 web search hook에 근거한 세 가지 personalized format(email, call script, LinkedIn message)을 제공합니다. 가능하면 사용자의 company positioning에 맞춥니다.

**Prospecting** — Net-new company(ProspectorOrganization)와 existing account(Organization)를 구분합니다. Iterative refinement와 "find companies like [X]" 같은 lookalike search를 지원합니다. Web search가 net-new result를 enrich합니다.

**Weekly Prep Brief** — 다음 7일의 모든 external call을 다루는 full briefing입니다. Company snapshot, attendee profile, signal, meeting별 recommended objective를 포함합니다.

## 설정

1. Cowork settings에서 Common Room MCP server가 connected 및 authenticated 상태인지 확인합니다.
2. 선택 사항으로 call prep 및 weekly briefing의 자동 meeting lookup을 위해 calendar MCP server를 연결합니다.
3. 이 plugin을 설치합니다. 모든 skill과 command를 즉시 사용할 수 있습니다.

## User context

User territory로 scope하는 모든 skill은 Common Room에서 `Me` object를 자동으로 가져옵니다. 이를 통해 user profile, role, "My Segments"를 얻고 query가 기본적으로 해당 territory를 사용하게 합니다. 자세한 내용은 `references/me-context.md`를 참고하세요.

Company context가 있으면 skill이 user product와 ICP에 맞게 recommendation을 조정합니다. 자세한 내용은 `references/my-company-context.md`를 참고하세요.

## 커스터마이징

Calendar connector와 tool reference 동작 방식은 `CONNECTORS.md`를 참고하세요.
