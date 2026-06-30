# Engineering 플러그인

Software engineering plugin입니다. 주로 Anthropic의 agentic desktop application인 [Cowork](https://claude.com/product/cowork)용으로 설계되었지만 Claude Code에서도 동작합니다. Standup, code review, architecture decision, incident response, debugging, technical documentation을 돕습니다. 어떤 engineering team에서도 사용할 수 있으며, 입력만으로도 standalone 동작하고 source control, project tracker, monitoring tool을 연결하면 더 강력해집니다.

## 설치

```bash
claude plugins add knowledge-work-plugins/engineering
```

## 명령

Slash command로 호출하는 명시적 workflow입니다:

| 명령 | 설명 |
|---|---|
| `/standup` | Commit, PR, ticket, chat 등 recent activity에서 standup update를 생성합니다. |
| `/review` | Code change를 security, performance, style, correctness 관점에서 review합니다. |
| `/debug` | Reproduce, isolate, diagnose, fix로 이어지는 structured debugging session을 실행합니다. |
| `/architecture` | Trade-off analysis가 포함된 ADR format으로 architecture decision을 작성하거나 평가합니다. |
| `/incident` | Triage, communication, mitigation, postmortem 작성을 포함한 incident response workflow를 실행합니다. |
| `/deploy-checklist` | Test 검증, change review, dependency check, rollback plan 확인을 포함한 pre-deployment checklist를 생성합니다. |

모든 명령은 code paste, system description, file upload만으로 **standalone** 동작하며, MCP connector를 연결하면 더 강력해집니다.

## 스킬

관련 상황에서 Claude가 자동으로 사용하는 domain knowledge입니다:

| 스킬 | 설명 |
|---|---|
| `code-review` | 코드 변경을 security, performance, correctness 관점에서 review합니다. PR URL 또는 diff, "review this before I merge", "is this code safe?" 또는 N+1 query, injection risk, missing edge case, error handling gap을 확인할 때 트리거됩니다. |
| `incident-response` | Incident response 워크플로를 실행합니다. Triage, communication, postmortem 작성을 다룹니다. "we have an incident", "production is down", severity assessment가 필요한 alert, incident 중 status update, 해결 후 blameless postmortem 작성 시 트리거됩니다. |
| `system-design` | 시스템, 서비스, architecture를 설계합니다. "design a system for", "how should we architect", "system design for", "what's the right architecture for" 또는 API design, data modeling, service boundary에 도움이 필요할 때 트리거됩니다. |
| `tech-debt` | Technical debt를 식별, 분류, 우선순위화합니다. "tech debt", "technical debt audit", "what should we refactor", "code health" 또는 code quality, refactoring priority, maintenance backlog에 대한 질문에서 트리거됩니다. |
| `testing-strategy` | 테스트 전략과 테스트 계획을 설계합니다. "how should we test", "test strategy for", "write tests for", "test plan", "what tests do we need" 또는 testing approach, coverage, test architecture에 도움이 필요할 때 트리거됩니다. |
| `documentation` | 기술 문서를 작성하고 유지합니다. "write docs for", "document this", "create a README", "write a runbook", "onboarding guide" 또는 API docs, architecture docs, operational runbook 등 기술 글쓰기에 도움이 필요할 때 트리거됩니다. |

## 예시 워크플로

### Morning Standup

```
/standup
```

If your tools are connected, I'll pull your recent commits, PR activity, and ticket updates. Otherwise, tell me what you worked on and I'll format it.

### Code Review

```
/review https://github.com/org/repo/pull/123
```

Share a PR link, paste a diff, or point to files. Get a structured review covering security, performance, correctness, and style.

### Debugging an Issue

```
/debug Users are getting 500 errors on the checkout page
```

Walk through a structured debugging process: reproduce, isolate, diagnose, fix. I'll help you think through it systematically.

### Architecture Decision

```
/architecture Should we use a message queue or direct API calls between services?
```

Get a structured ADR with options analysis, trade-offs, and a recommendation.

### Incident Response

```
/incident The payments service is returning 503s
```

Start an incident workflow: triage severity, draft communications, track timeline, and generate a postmortem when resolved.

### Pre-Deploy Check

```
/deploy-checklist auth-service v2.3.0
```

Get a customized deployment checklist based on your service and what's changing.

## 단독 사용과 통합 사용

모든 명령과 스킬은 integration 없이도 동작합니다:

| 할 수 있는 일 | 단독 사용 | 통합 사용 시 |
|-----------------|------------|-------------------|
| Standup updates | Describe your work | Source control, Project tracker, Chat |
| Code review | Paste diff or code | Source control (pull PRs automatically) |
| Debug sessions | Describe the problem | Monitoring (pull logs and metrics) |
| Architecture decisions | Describe the system | Knowledge base (find prior ADRs) |
| Incident response | Describe the incident | Monitoring, Incident management, Chat |
| Deploy checklists | Describe the deploy | CI/CD, Source control |

## MCP 통합

> If you see unfamiliar placeholders or need to check which tools are connected, see [CONNECTORS.md](CONNECTORS.md).

더 풍부한 경험을 위해 tool을 연결하세요:

| 범주 | 예시 | 가능해지는 일 |
|---|---|---|
| **Source control** | GitHub, GitLab | PR diffs, commit history, branch status |
| **Project tracker** | Linear, Jira, Asana | Ticket status, sprint data, assignments |
| **Monitoring** | Datadog, New Relic | Logs, metrics, alerts, dashboards |
| **Incident management** | PagerDuty, Opsgenie | On-call schedules, incident tracking, paging |
| **Chat** | Slack, Teams | Team discussions, standup channels |
| **Knowledge base** | Notion, Confluence | ADRs, runbooks, onboarding docs |

지원되는 integration 전체 목록은 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

## 설정

개인화를 위해 `engineering/.claude/settings.local.json`에 local settings file을 만드세요:

```json
{
  "name": "Your Name",
  "title": "Software Engineer",
  "team": "Your Team",
  "company": "Your Company",
  "techStack": ["Python", "TypeScript", "PostgreSQL", "AWS"],
  "defaultBranch": "main",
  "deployProcess": "canary"
}
```

설정되어 있지 않으면 플러그인이 이 정보를 대화형으로 물어봅니다.
