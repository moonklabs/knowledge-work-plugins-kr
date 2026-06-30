# Common Room 플러그인

Common Room 기반 GTM workflow입니다. Account research, contact research, call prep, personalized outreach, prospecting, weekly briefing을 지원합니다.

## Overview

이 plugin은 Claude를 Common Room MCP server에 연결하고, 가장 일반적인 rep workflow를 다루는 여섯 가지 skill을 제공합니다. 모든 output은 실제 Common Room signal data에 근거합니다. 1st-party product signal, 2nd-party community signal, 3rd-party intent signal, RoomieAI와 Spark의 enrichment를 활용합니다.

## Requirements

- **Common Room MCP** (`mcp.commonroom.io/mcp`) must be connected and authenticated. This is the primary data source for all plugin functionality.
- **Calendar connector** (optional) — enables automatic meeting lookup in `call-prep` and `weekly-prep-brief`. If not connected, both skills ask the user for meeting details instead.

## 스킬

Skills are triggered conversationally. Describe what you want and Claude will load the right skill automatically.

| Skill | Trigger phrases |
|-------|----------------|
| `account-research` | Common Room data를 사용해 company를 조사합니다. 'research [company]', 'tell me about [domain]', 'pull up signals for [account]', 'what's going on with [company]' 또는 account-level question에서 트리거됩니다. |
| `contact-research` | Common Room data를 사용해 특정 person을 조사합니다. 'who is [name]', 'look up [email]', 'research [contact]', 'is [name] a warm lead' 또는 contact-level question에서 트리거됩니다. |
| `call-prep` | Common Room signal을 사용해 customer 또는 prospect call을 준비합니다. 'prep me for my call with [company]', 'prepare for a meeting with [company]', 'what should I know before talking to [company]' 또는 call preparation request에서 트리거됩니다. |
| `compose-outreach` | Common Room signal을 사용해 개인화된 outreach message를 생성합니다. 'draft outreach to [person]', 'write an email to [name]', 'compose a message for [contact]' 또는 outreach drafting request에서 트리거됩니다. |
| `prospect` | Common Room Prospector를 사용해 targeted account 또는 contact list를 만듭니다. 'find companies that match [criteria]', 'build a prospect list', 'find contacts at [type of company]', 'show me companies hiring [role]' 또는 list-building request에서 트리거됩니다. |
| `weekly-prep-brief` | 다음 7일 동안의 모든 external call에 대한 comprehensive weekly briefing을 생성합니다. 'weekly prep brief', 'prepare my week', 'what calls do I have this week', 'Monday prep' 또는 weekly planning request에서 트리거됩니다. |

## 명령

Two commands for complex workflows that benefit from explicit invocation:

| Command | Usage |
|---------|-------|
| `/generate-account-plan <company>` | Comprehensive strategic account plan with stakeholder mapping, engagement analysis, opportunities, risks, and action items |
| `/weekly-brief [date range]` | Generate a full weekly prep briefing (defaults to next 7 days) |

## What Each Skill Produces

**Account Research** — Handles four patterns: full overviews, targeted field questions, honest sparse-data responses, and combined MCP data + LLM reasoning. Includes web search for recent news. Automatically scopes to "My Segments."

**Contact Research** — Lookup by email, name+company, or social handle. Returns enriched identity, CRM fields, scores, website visits, activity history, Spark analyses, and conversation starters.

**Call Prep** — Company snapshot, per-attendee profiles, signal highlights, tailored talking points, likely objections, and recommended call outcome. Prioritizes Gong/call recording activity. Calendar-aware if connected.

**Compose Outreach** — Three personalized formats (email, call script, LinkedIn message) grounded in Common Room signals and web search hooks. Tailored to the user's company positioning when available.

**Prospecting** — Distinguishes between net-new companies (ProspectorOrganization) and existing accounts (Organization). Supports iterative refinement and lookalike search ("find companies like [X]"). Web search enriches net-new results.

**Weekly Prep Brief** — Full briefing covering every external call in the next 7 days: company snapshot, attendee profiles, signals, and recommended objectives per meeting.

## Setup

1. Ensure the Common Room MCP server is connected and authenticated in your Cowork settings.
2. (Optional) Connect a calendar MCP server for automatic meeting lookup in call prep and weekly briefings.
3. Install this plugin. All skills and commands are available immediately.

## User Context

All skills that scope to a user's territory automatically fetch the `Me` object from Common Room. This provides the user's profile, role, and "My Segments" — ensuring queries default to their territory. See `references/me-context.md` for details.

When company context is available, skills tailor recommendations to the user's product and ICP. See `references/my-company-context.md` for details.

## Customization

See `CONNECTORS.md` for details on the calendar connector and how tool references work.
