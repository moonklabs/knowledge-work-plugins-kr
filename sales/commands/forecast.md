---
description: 베스트/가능/최악 시나리오, 커밋 vs 업사이드, 갭 분석을 포함한 가중 세일즈 예측 생성
argument-hint: "<period>"
---

# /forecast

> 낯선 플레이스홀더가 보이거나 어떤 도구가 연결되어 있는지 확인하려면 [CONNECTORS.md](../CONNECTORS.md)를 보세요.

리스크 분석과 커밋 추천이 포함된 가중 세일즈 예측을 생성합니다.

## 사용법

```
/forecast
```

그 다음 파이프라인 데이터와 목표를 제공하세요.

---

## 동작 방식

```
┌─────────────────────────────────────────────────────────────────┐
│                        FORECAST                                  │
├─────────────────────────────────────────────────────────────────┤
│  STANDALONE (always works)                                       │
│  ✓ CRM에서 CSV 내보내기 업로드                                  │
│  ✓ 또는 파이프라인 딜을 붙여넣기/설명                            │
│  ✓ 쿼터와 타임라인 설정                                          │
│  ✓ 단계 확률 기반 가중 예측                                     │
│  ✓ 리스크 반영 전망(베스트/가능/최악)                           │
│  ✓ 커밋 vs 업사이드 분해                                        │
│  ✓ 갭 분석 및 추천                                               │
├─────────────────────────────────────────────────────────────────┤
│  SUPERCHARGED (when you connect your tools)                      │
│  + CRM: 파이프라인 자동 수집, 실시간 데이터                     │
│  + 단계/세그먼트/딜 규모별 과거 승률                             │
│  + 리스크 스코어링용 활동 신호                                   │
│  + 자동 갱신 및 추적                                             │
└─────────────────────────────────────────────────────────────────┘
```

---

## 필요한 정보

### 1단계: 파이프라인 데이터

**옵션 A: CSV 업로드**
CRM(예: Salesforce, HubSpot)에서 파이프라인을 내보내세요. 최소한 다음이 필요합니다.
- Deal/Opportunity name
- Amount
- Stage
- Close date

있으면 도움이 되는 항목:
- Owner (if team forecast)
- Last activity date
- Created date
- Account name

**옵션 B: 딜 붙여넣기**
```
Acme Corp - $50K - Negotiation - closes Jan 31
TechStart - $25K - Demo scheduled - closes Feb 15
BigCo - $100K - Discovery - closes Mar 30
```

**옵션 C: 테리토리 설명**
"I have 8 deals in pipeline totaling $400K. Two are in negotiation ($120K), three in evaluation ($180K), three in discovery ($100K)."

### 2단계: 목표

- **Quota**: 목표 수치는 무엇인가요? (예: "이번 분기 $500K")
- **Timeline**: 기간 종료일은 언제인가요? (예: "Q1은 3월 31일 종료")
- **Already closed**: 이번 기간에 이미 확정한 금액은 얼마인가요?

---

## 출력

```markdown
# Sales Forecast: [Period]

**Generated:** [Date]
**Data Source:** [CSV upload / Manual input / CRM]

---

## Summary

| Metric | Value |
|--------|-------|
| **Quota** | $[X] |
| **Closed to Date** | $[X] ([X]% of quota) |
| **Open Pipeline** | $[X] |
| **Weighted Forecast** | $[X] |
| **Gap to Quota** | $[X] |
| **Coverage Ratio** | [X]x |

---

## Forecast Scenarios

| Scenario | Amount | % of Quota | Assumptions |
|----------|--------|------------|-------------|
| **Best Case** | $[X] | [X]% | All deals close as expected |
| **Likely Case** | $[X] | [X]% | Stage-weighted probabilities |
| **Worst Case** | $[X] | [X]% | Only commit deals close |

---

## Pipeline by Stage

| Stage | # Deals | Total Value | Probability | Weighted Value |
|-------|---------|-------------|-------------|----------------|
| Negotiation | [X] | $[X] | 80% | $[X] |
| Proposal | [X] | $[X] | 60% | $[X] |
| Evaluation | [X] | $[X] | 40% | $[X] |
| Discovery | [X] | $[X] | 20% | $[X] |
| **Total** | [X] | $[X] | — | $[X] |

---

## Commit vs. Upside

### Commit (High Confidence)
확신할 수 있는 딜:

| Deal | Amount | Stage | Close Date | Why Commit |
|------|--------|-------|------------|------------|
| [Deal] | $[X] | [Stage] | [Date] | [Reason] |

**Total Commit:** $[X]

### Upside (Lower Confidence)
닫힐 수 있지만 리스크가 있는 딜:

| Deal | Amount | Stage | Close Date | Risk Factor |
|------|--------|-------|------------|-------------|
| [Deal] | $[X] | [Stage] | [Date] | [Risk] |

**Total Upside:** $[X]

---

## Risk Flags

| Deal | Amount | Risk | Recommendation |
|------|--------|------|----------------|
| [Deal] | $[X] | Close date passed | Update close date or move to lost |
| [Deal] | $[X] | No activity in 14+ days | Re-engage or downgrade stage |
| [Deal] | $[X] | Close date this week, still in discovery | Unlikely to close — push out |

---

## Gap Analysis

**To hit quota, you need:** $[X] more

**Options to close the gap:**
1. **Accelerate [Deal]** — Currently [stage], worth $[X]. If you can close by [date], you're at [X]% of quota.
2. **Revive [Stalled Deal]** — Last active [date]. Worth $[X]. Reach out to [contact].
3. **New pipeline needed** — You need $[X] in new opportunities at [X]x coverage to be safe.

---

## Recommendations

1. [ ] [Specific action for highest-impact deal]
2. [ ] [Action for at-risk deal]
3. [ ] [Pipeline generation recommendation if gap exists]
```

---

## Stage Probabilities (Default)

If you don't provide custom probabilities, I'll use:

| Stage | Default Probability |
|-------|---------------------|
| Closed Won | 100% |
| Negotiation / Contract | 80% |
| Proposal / Quote | 60% |
| Evaluation / Demo | 40% |
| Discovery / Qualification | 20% |
| Prospecting / Lead | 10% |

단계나 확률이 다르다면 알려 주세요.

---

## CRM 연결 시

- 파이프라인을 자동으로 가져옵니다
- 실제 과거 승률을 사용합니다
- 활동 최신성을 리스크 점수에 반영합니다
- Track forecast changes over time
- Compare to previous forecasts

---

## Tips

1. **Be honest about commit** — Only commit deals you'd bet on. Upside is for everything else.
2. **Update close dates** — Stale close dates kill forecast accuracy. Push out deals that won't close in time.
3. **Coverage matters** — 3x pipeline coverage is healthy. Below 2x is risky.
4. **Activity = signal** — Deals with no recent activity are at higher risk than stage suggests.
