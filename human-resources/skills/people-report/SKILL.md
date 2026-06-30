---
name: people-report
description: 인원수, 이탈률, 다양성, 조직 건강 보고서를 생성합니다. 리더십용 인원수 스냅샷, team별 이직 추세 분석, 다양성 대표성 지표 준비, 조직 전반의 관리 범위과 이탈 위험 평가에 사용합니다.
argument-hint: "<report type — headcount, attrition, diversity, org health 등>"
---

# /people-report

> If you see unfamiliar placeholders or need to check which tools are connected, see [CONNECTORS.md](../../CONNECTORS.md).

Generate people analytics reports from your HR data. Analyze workforce data to surface trends, risks, and opportunities.

## Usage

```
/people-report $ARGUMENTS
```

## Report Types

**Headcount**: Current org snapshot — by team, location, level, tenure
**Attrition**: Turnover analysis — voluntary/involuntary, by team, trends
**Diversity**: Representation metrics — by level, team, pipeline
**Org Health**: Span of control, management layers, team sizes, flight risk

## Key Metrics

### Retention
- Overall attrition rate (voluntary + involuntary)
- Regrettable attrition rate
- Average tenure
- Flight risk indicators

### Diversity
- Representation by level, team, and function
- Pipeline diversity (hiring funnel by demographic)
- Promotion rates by group
- Pay equity analysis

### Engagement
- Survey scores and trends
- eNPS (Employee Net Promoter Score)
- Participation rates
- Open-ended feedback themes

### Productivity
- Revenue per employee
- Span of control efficiency
- Time to productivity for new hires

## Approach

1. Understand what question they're trying to answer
2. Identify the right data (upload, paste, or pull from ~~HRIS)
3. Analyze with appropriate statistical methods
4. Present findings with context and caveats
5. Recommend specific actions based on data

## What I Need From You

Upload a CSV or describe your data. Helpful fields:
- Employee name/ID, department, team
- Title, level, location
- Start date, end date (if applicable)
- Manager, compensation (if relevant)
- Demographics (for diversity reports, if available)

## Output

```markdown
## People Report: [Type] — [Date]

### Executive Summary
[2-3 key takeaways]

### Key Metrics
| Metric | Value | Trend |
|--------|-------|-------|
| [Metric] | [Value] | [up/down/flat] |

### Detailed Analysis
[Charts, tables, and narrative for the specific report type]

### Recommendations
- [Data-driven recommendation]
- [Action item]

### Methodology
[How the numbers were calculated, any caveats]
```

## If Connectors Available

If **~~HRIS** is connected:
- Pull live employee data — headcount, tenure, department, level
- Generate reports without needing a CSV upload

If **~~chat** is connected:
- Offer to share the report summary in a relevant channel
