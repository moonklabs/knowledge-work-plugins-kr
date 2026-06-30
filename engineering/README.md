# Engineering 플러그인

Software engineering 플러그인입니다. 주로 Anthropic의 agentic 데스크톱 애플리케이션인 [Cowork](https://claude.com/product/cowork)용으로 설계되었지만 Claude Code에서도 동작합니다. Standup, code review, architecture decision, incident response, debugging, technical documentation을 돕습니다. 어떤 engineering team에서도 사용할 수 있으며, 입력만으로도 단독 동작하고 source control, project tracker, monitoring tool을 연결하면 더 강력해집니다.

## 설치

```bash
claude plugins add knowledge-work-plugins/engineering
```

## 명령

Slash command로 호출하는 명시적 워크플로입니다:

| 명령 | 설명 |
|---|---|
| `/standup` | Commit, PR, ticket, chat 등 recent activity에서 standup update를 생성합니다. |
| `/review` | Code change를 security, performance, style, correctness 관점에서 검토합니다. |
| `/debug` | 재현, 격리, 진단, 수정으로 이어지는 구조화된 debugging session을 실행합니다. |
| `/architecture` | Trade-off analysis가 포함된 ADR format으로 architecture decision을 작성하거나 평가합니다. |
| `/incident` | Triage, communication, mitigation, postmortem 작성을 포함한 incident response 워크플로를 실행합니다. |
| `/deploy-checklist` | Test 검증, change review, dependency check, rollback plan 확인을 포함한 pre-deployment checklist를 생성합니다. |

모든 명령은 code paste, system description, file upload만으로 **단독** 동작하며, MCP connector를 연결하면 더 강력해집니다.

## 스킬

관련 상황에서 Claude가 자동으로 사용하는 domain knowledge입니다:

| 스킬 | 설명 |
|---|---|
| `code-review` | 코드 변경을 security, performance, correctness 관점에서 검토합니다. PR URL 또는 diff, "review this before I merge", "is this code safe?" 또는 N+1 query, injection risk, missing edge case, error handling gap을 확인할 때 트리거됩니다. |
| `incident-response` | Incident response 워크플로를 실행합니다. Triage, communication, postmortem 작성을 다룹니다. "we have an incident", "production is down", severity assessment가 필요한 alert, incident 중 status update, 해결 후 blameless postmortem 작성 시 트리거됩니다. |
| `system-design` | 시스템, 서비스, architecture를 설계합니다. "design a system for", "how should we architect", "system design for", "what's the right architecture for" 또는 API design, data modeling, service boundary에 도움이 필요할 때 트리거됩니다. |
| `tech-debt` | Technical debt를 식별, 분류, 우선순위화합니다. "tech debt", "technical debt audit", "what should we refactor", "code health" 또는 code quality, refactoring priority, maintenance backlog에 대한 질문에서 트리거됩니다. |
| `testing-strategy` | 테스트 전략과 테스트 계획을 설계합니다. "how should we test", "test strategy for", "write tests for", "test plan", "what tests do we need" 또는 testing approach, coverage, test architecture에 도움이 필요할 때 트리거됩니다. |
| `documentation` | 기술 문서를 작성하고 유지합니다. "write docs for", "document this", "create a README", "write a runbook", "onboarding guide" 또는 API docs, architecture docs, operational runbook 등 기술 글쓰기에 도움이 필요할 때 트리거됩니다. |

## 예시 워크플로

### Morning standup

```text
/standup
```

Tool이 연결되어 있으면 recent commit, PR activity, ticket update를 가져옵니다. 그렇지 않으면 작업한 내용을 알려주면 standup format으로 정리합니다.

### Code review

```text
/review https://github.com/org/repo/pull/123
```

PR link를 공유하거나 diff를 붙여넣거나 file을 지정하세요. Security, performance, correctness, style을 다루는 structured review를 받습니다.

### Issue debugging

```text
/debug 사용자가 checkout page에서 500 error를 보고 있습니다
```

Reproduce, isolate, diagnose, fix로 이어지는 structured debugging process를 진행합니다. Systematic하게 생각하도록 돕습니다.

### Architecture decision

```text
/architecture 서비스 간 통신에 message queue를 써야 할까요, direct API call을 써야 할까요?
```

Option analysis, trade-off, recommendation이 포함된 structured ADR을 받습니다.

### Incident response

```text
/incident payments service가 503을 반환하고 있습니다
```

Incident workflow를 시작합니다. Severity를 triage하고 communication 초안을 작성하며 timeline을 track하고 resolved 후 postmortem을 생성합니다.

### Pre-deploy check

```text
/deploy-checklist auth-service v2.3.0
```

Service와 변경 내용을 기준으로 customized deployment checklist를 받습니다.

## 단독 사용과 통합 사용

모든 명령과 스킬은 integration 없이도 동작합니다:

| 할 수 있는 일 | 단독 사용 | 통합 사용 시 |
|-----------------|------------|-------------------|
| Standup update | 작업 내용 설명 | Source control, Project tracker, Chat |
| Code review | Diff 또는 code 붙여넣기 | Source control(PR 자동 수집) |
| Debug session | Problem 설명 | Monitoring(log와 metric 수집) |
| Architecture decision | System 설명 | Knowledge base(prior ADR 검색) |
| Incident response | Incident 설명 | Monitoring, Incident management, Chat |
| Deploy checklist | Deploy 설명 | CI/CD, Source control |

## MCP 통합

> 익숙하지 않은 placeholder가 보이거나 연결된 도구를 확인해야 한다면 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

더 풍부한 경험을 위해 tool을 연결하세요:

| 범주 | 예시 | 가능해지는 일 |
|---|---|---|
| **Source control** | GitHub, GitLab | PR diff, commit history, branch status |
| **Project tracker** | Linear, Jira, Asana | Ticket status, sprint data, assignment |
| **Monitoring** | Datadog, New Relic | Log, metric, alert, dashboard |
| **Incident management** | PagerDuty, Opsgenie | On-call schedule, incident tracking, paging |
| **Chat** | Slack, Teams | Team discussion, standup channel |
| **Knowledge base** | Notion, Confluence | ADR, runbook, onboarding doc |

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
