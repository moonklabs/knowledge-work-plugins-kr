---
description: 파이프라인 상태 분석, 딜 우선순위 지정, 리스크 플래그, 주간 액션 플랜 생성
argument-hint: "<segment or rep>"
---

# /pipeline-review

> 낯선 플레이스홀더가 보이거나 어떤 도구가 연결되어 있는지 확인하려면 [CONNECTORS.md](../CONNECTORS.md)를 보세요.

파이프라인 건전성을 분석하고 딜 우선순위를 정해 어디에 집중해야 할지 실행 가능한 추천을 제공합니다.

## 사용법

```
/pipeline-review
```

그 다음 파이프라인 데이터를 제공하세요.

---

## 동작 방식

```
┌─────────────────────────────────────────────────────────────────┐
│                     PIPELINE REVIEW                              │
├─────────────────────────────────────────────────────────────────┤
│  STANDALONE (always works)                                       │
│  ✓ CRM CSV 내보내기 업로드                                      │
│  ✓ 또는 딜을 붙여넣거나 설명                                    │
│  ✓ 상태 점검: 오래됨/정체/고위험 딜 플래그                     │
│  ✓ 우선순위: 영향도와 클로즈 가능성 기반 랭킹                  │
│  ✓ 위생 점검: 데이터 누락, 잘못된 종료일, 단일 접점            │
│  ✓ 주간 액션 플랜: 이번 주 집중 항목 제시                      │
├─────────────────────────────────────────────────────────────────┤
│  SUPERCHARGED (when you connect your tools)                      │
│  + CRM: 파이프라인 자동 수집 및 레코드 업데이트                 │
│  + 활동 데이터 기반 참여도 스코어링                             │
│  + 과거 패턴 기반 리스크 예측                                   │
│  + Calendar: 딜별 예정 미팅 확인                                │
└─────────────────────────────────────────────────────────────────┘
```

---

## 필요한 정보

**옵션 A: CSV 업로드**
CRM(예: Salesforce, HubSpot)에서 파이프라인을 내보내세요. 다음 필드가 있으면 좋습니다.
- Deal/Opportunity name
- Account name
- Amount
- Stage
- Close date
- Created date
- Last activity date
- Owner (팀 리뷰 시)
- Primary contact

**옵션 B: 딜 붙여넣기**
```
Acme Corp - $50K - Negotiation - closes Jan 31 - last activity Jan 20
TechStart - $25K - Demo scheduled - closes Feb 15 - no activity in 3 weeks
BigCo - $100K - Discovery - closes Mar 30 - created last week
```

**옵션 C: 파이프라인 설명**
"I have 12 deals. Two big ones in negotiation that I'm confident about. Three stuck in discovery for over a month. The rest are mid-stage but I haven't talked to some of them in a while."

---

## 출력

```markdown
# Pipeline Review: [Date]

**Data Source:** [CSV upload / Manual input / CRM]
**Deals Analyzed:** [X]
**Total Pipeline Value:** $[X]

---

## Pipeline Health Score: [X/100]

| Dimension | Score | Issue |
|-----------|-------|-------|
| **Stage Progression** | [X]/25 | [X] deals stuck in same stage 30+ days |
| **Activity Recency** | [X]/25 | [X] deals with no activity in 14+ days |
| **Close Date Accuracy** | [X]/25 | [X] deals with close date in past |
| **Contact Coverage** | [X]/25 | [X] deals single-threaded |

---

## Priority Actions This Week

### 1. [Highest Priority Deal]
**Why:** [Reason — large, closing soon, at risk, etc.]
**Action:** [Specific next step]
**Impact:** $[X] if you close it

### 2. [Second Priority]
**Why:** [Reason]
**Action:** [Next step]

### 3. [Third Priority]
**Why:** [Reason]
**Action:** [Next step]

---

## Deal Prioritization Matrix

### Close This Week (Focus Time Here)
| Deal | Amount | Stage | Close Date | Next Action |
|------|--------|-------|------------|-------------|
| [Deal] | $[X] | [Stage] | [Date] | [Action] |

### Close This Month (Keep Warm)
| Deal | Amount | Stage | Close Date | Status |
|------|--------|-------|------------|--------|
| [Deal] | $[X] | [Stage] | [Date] | [Status] |

### Nurture (Check-in Periodically)
| Deal | Amount | Stage | Close Date | Status |
|------|--------|-------|------------|--------|
| [Deal] | $[X] | [Stage] | [Date] | [Status] |

---

## Risk Flags

### Stale Deals (No Activity 14+ Days)
| Deal | Amount | Last Activity | Days Silent | Recommendation |
|------|--------|---------------|-------------|----------------|
| [Deal] | $[X] | [Date] | [X] | [Re-engage / Downgrade / Remove] |

### Stuck Deals (Same Stage 30+ Days)
| Deal | Amount | Stage | Days in Stage | Recommendation |
|------|--------|-------|---------------|----------------|
| [Deal] | $[X] | [Stage] | [X] | [Push / Multi-thread / Qualify out] |

### Past Close Date
| Deal | Amount | Close Date | Days Overdue | Recommendation |
|------|--------|------------|--------------|----------------|
| [Deal] | $[X] | [Date] | [X] | [Update date / Push to next quarter / Close lost] |

### Single-Threaded (Only One Contact)
| Deal | Amount | Contact | Risk | Recommendation |
|------|--------|---------|------|----------------|
| [Deal] | $[X] | [Name] | Champion leaves = deal dies | [Identify additional stakeholders] |

---

## Hygiene Issues

| Issue | Count | Deals | Action |
|-------|-------|-------|--------|
| Missing close date | [X] | [List] | Add realistic close dates |
| Missing amount | [X] | [List] | Estimate or qualify |
| Missing next step | [X] | [List] | Define next action |
| No primary contact | [X] | [List] | Assign contact |

---

## Pipeline Shape

### By Stage
| Stage | # Deals | Value | % of Pipeline |
|-------|---------|-------|---------------|
| [Stage] | [X] | $[X] | [X]% |

### By Close Month
| Month | # Deals | Value |
|-------|---------|-------|
| [Month] | [X] | $[X] |

### By Deal Size
| Size | # Deals | Value |
|------|---------|-------|
| $100K+ | [X] | $[X] |
| $50K-100K | [X] | $[X] |
| $25K-50K | [X] | $[X] |
| <$25K | [X] | $[X] |

---

## Recommendations

### This Week
1. [ ] [Specific action for priority deal 1]
2. [ ] [Action for at-risk deal]
3. [ ] [Hygiene task]

### This Month
1. [ ] [Strategic action]
2. [ ] [Pipeline building if needed]

---

## Deals to Consider Removing

These deals may be dead weight:

| Deal | Amount | Reason | Recommendation |
|------|--------|--------|----------------|
| [Deal] | $[X] | [No activity 60+ days, no response] | Mark closed-lost |
| [Deal] | $[X] | [Pushed 3+ times, no champion] | Qualify out |
```

---

## 우선순위 프레임워크

다음 기준으로 딜을 랭킹합니다.

| 요소 | 가중치 | 기준 |
|--------|--------|-----------------|
| **Close Date** | 30% | 마감이 임박한 딜 우선 |
| **Deal Size** | 25% | 금액이 큰 딜 우선 |
| **Stage** | 20% | 후반 단계일수록 우선 |
| **Activity** | 15% | 최근 활동 있는 딜 우선 |
| **Risk** | 10% | 리스크 낮은 딜 우선 |

가중치는 조정할 수 있습니다. 예: "빠른 성과가 필요하니 Close Date를 더 높게".

---

## CRM 연결 시

- 파이프라인을 자동으로 수집합니다
- 종료일, 단계, 다음 단계를 레코드에 업데이트합니다
- 후속 작업을 생성합니다
- 데이터 위생 개선 추이를 추적합니다

---

## 팁

1. **주간 리뷰**: 파이프라인 건전성은 빠르게 악화됩니다. 주간 점검이 조기 대응에 효과적입니다.
2. **죽은 딜 정리**: 오래된 기회는 파이프라인을 부풀리고 예측을 왜곡합니다.
3. **멀티 스레딩**: 핵심 연락처 1명에 의존하지 말고 복수 이해관계자를 확보하세요.
4. **종료일 신뢰성**: 종료일은 "희망일"이 아니라 실제 서명 예상일이어야 합니다.
