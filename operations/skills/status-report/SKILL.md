---
name: status-report
description: KPI, risk, action item이 포함된 status report를 생성합니다. Leadership용 weekly/monthly update 작성, green/yellow/red status로 project health 요약, stakeholder attention이 필요한 risk와 decision surface, project tracker activity를 readable narrative로 변환할 때 사용합니다.
argument-hint: "[weekly | monthly | quarterly] [project or team]"
---

# /status-report

> If you see unfamiliar placeholders or need to check which tools are connected, see [CONNECTORS.md](../../CONNECTORS.md).

Generate a polished status report for leadership or stakeholders. See the **risk-assessment** skill for risk matrix frameworks and severity definitions.

## Usage

```
/status-report $ARGUMENTS
```

## Output

```markdown
## Status Report: [Project/Team] — [Period]
**Author:** [Name] | **Date:** [Date]

### Executive Summary
[3-4 sentence overview — what's on track, what needs attention, key wins]

### Overall Status: 🟢 On Track / 🟡 At Risk / 🔴 Off Track

### Key Metrics
| Metric | Target | Actual | Trend | Status |
|--------|--------|--------|-------|--------|
| [KPI] | [Target] | [Actual] | [up/down/flat] | 🟢/🟡/🔴 |

### Accomplishments This Period
- [Win 1]
- [Win 2]

### In Progress
| Item | Owner | Status | ETA | Notes |
|------|-------|--------|-----|-------|
| [Item] | [Person] | [Status] | [Date] | [Context] |

### Risks and Issues
| Risk/Issue | Impact | Mitigation | Owner |
|------------|--------|------------|-------|
| [Risk] | [Impact] | [What we're doing] | [Who] |

### Decisions Needed
| Decision | Context | Deadline | Recommended Action |
|----------|---------|----------|--------------------|
| [Decision] | [Why it matters] | [When] | [What I recommend] |

### Next Period Priorities
1. [Priority 1]
2. [Priority 2]
3. [Priority 3]
```

## If Connectors Available

If **~~project tracker** is connected:
- Pull project status, completed items, and upcoming milestones automatically
- Identify at-risk items and overdue tasks

If **~~chat** is connected:
- Scan recent team discussions for decisions and blockers to include
- Offer to post the finished report to a channel

If **~~calendar** is connected:
- Reference key meetings and decisions from the reporting period

## Tips

1. **Lead with the headline** — Busy leaders read the first 3 lines. Make them count.
2. **Be honest about risks** — Surfacing issues early builds trust. Surprises erode it.
3. **Make decisions easy** — For each decision needed, provide context and a recommendation.
