# Operations 플러그인

Business operations plugin입니다. 주로 Anthropic의 agentic desktop application인 [Cowork](https://claude.com/product/cowork)용으로 설계되었지만 Claude Code에서도 동작합니다. Vendor management, process documentation, change management, capacity planning, compliance tracking, resource planning을 돕습니다. 어떤 ops team에서도 사용할 수 있으며, 입력만으로 standalone 동작하고 ITSM, project tracker 및 다른 tool을 연결하면 더 강력해집니다.

## 설치

```bash
claude plugins add knowledge-work-plugins/operations
```

## 명령

Slash command로 호출하는 명시적 workflow입니다:

| 명령 | 설명 |
|---|---|
| `/vendor-review` | Cost analysis, risk assessment, contract summary, renewal recommendation으로 vendor를 평가합니다. |
| `/process-doc` | Flowchart, RACI matrix, SOP, runbook으로 business process를 문서화합니다. |
| `/change-request` | Impact analysis, rollback plan, approval routing이 포함된 change management request를 작성합니다. |
| `/capacity-plan` | Workload analysis, headcount modeling, utilization forecasting으로 resource capacity를 계획합니다. |
| `/status-report` | Leadership용 project update, KPI, risk, action item이 포함된 status report를 생성합니다. |
| `/runbook` | Recurring task를 위한 step-by-step operational runbook을 생성하거나 update합니다. |

All commands work **standalone** (provide context and details) and get **supercharged** with MCP connectors.

## 스킬

관련 상황에서 Claude가 자동으로 사용하는 domain knowledge입니다:

| 스킬 | 설명 |
|---|---|
| `vendor-management` | Evaluate, compare, and manage vendor relationships — contracts, performance, risk |
| `process-optimization` | Business process를 분석하고 개선합니다. "this process is slow", "how can we improve", "streamline this workflow", "too many steps", "bottleneck" 또는 사용자가 고치고 싶은 inefficient process를 설명할 때 트리거됩니다. |
| `change-management` | Plan and execute organizational or technical changes — communication, training, adoption |
| `risk-assessment` | Operational risk를 식별, 평가, 완화합니다. "what are the risks", "risk assessment", "risk register", "what could go wrong" 또는 project, vendor, process, decision 관련 risk를 평가할 때 트리거됩니다. |
| `compliance-tracking` | Compliance requirement와 audit readiness를 track합니다. "compliance", "audit prep", "SOC 2", "ISO 27001", "GDPR", "regulatory requirement" 또는 compliance activity tracking/preparation/documentation에 도움이 필요할 때 트리거됩니다. |
| `resource-planning` | Plan and optimize resource allocation — capacity, utilization, forecasting, budget |

## 예시 워크플로

### Evaluating a Vendor

```
/vendor-review
```

Provide the vendor name, contract details, or upload a proposal. Get a structured evaluation with cost analysis, risk flags, and a recommendation.

### Documenting a Process

```
/process-doc employee offboarding
```

Describe the process or walk me through it. Get a complete SOP with flowchart, RACI matrix, and step-by-step procedures.

### Submitting a Change Request

```
/change-request
```

Describe the change. Get an impact analysis, risk assessment, rollback plan, and communication template ready for approval.

### Planning Capacity

```
/capacity-plan
```

Upload team data or describe your resources. Get utilization analysis, bottleneck identification, and headcount recommendations.

### Leadership Status Report

```
/status-report
```

I'll pull updates from your connected tools (or ask you for input) and generate a polished status report with KPIs, risks, and next steps.

### Creating a Runbook

```
/runbook monthly close process
```

Walk me through the process once. I'll document it as a repeatable runbook with checklists, troubleshooting, and escalation paths.

## 단독 사용과 통합 사용

모든 명령과 스킬은 integration 없이도 동작합니다:

| 할 수 있는 일 | 단독 사용 | 통합 사용 시 |
|-----------------|------------|-------------------|
| Vendor reviews | Provide details, upload proposals | Procurement, Knowledge base |
| Process documentation | Describe the process | Knowledge base (existing docs) |
| Change requests | Describe the change | ITSM, Project tracker |
| Capacity planning | Upload data, describe team | Project tracker (workload data) |
| Status reports | Provide updates manually | Project tracker, Chat, Calendar |
| Runbooks | Walk through the process | Knowledge base, ITSM |

## MCP 통합

> If you see unfamiliar placeholders or need to check which tools are connected, see [CONNECTORS.md](CONNECTORS.md).

더 풍부한 경험을 위해 tool을 연결하세요:

| 범주 | 예시 | 가능해지는 일 |
|---|---|---|
| **ITSM** | ServiceNow, Zendesk | Ticket management, change requests, incident tracking |
| **Project tracker** | Asana, Jira, monday.com | Project status, resource allocation, task tracking |
| **Knowledge base** | Notion, Confluence | Process docs, runbooks, policies |
| **Chat** | Slack, Teams | Team coordination, approvals, status updates |
| **Calendar** | Google Calendar, Microsoft 365 | Meeting scheduling, deadline tracking |
| **Email** | Gmail, Microsoft 365 | Vendor communications, approvals |

지원되는 integration 전체 목록은 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

## 설정

Create a local settings file at `operations/.claude/settings.local.json` to personalize:

```json
{
  "company": "Your Company",
  "team": "Operations",
  "reportingCadence": "weekly",
  "approvalChain": ["Manager", "Director", "VP"],
  "complianceFrameworks": ["SOC 2", "ISO 27001"],
  "fiscalYearStart": "January"
}
```

설정되어 있지 않으면 플러그인이 이 정보를 대화형으로 물어봅니다.
