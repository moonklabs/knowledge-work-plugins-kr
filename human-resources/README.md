# HR 플러그인

People operations 플러그인입니다. 주로 Anthropic의 agentic 데스크톱 애플리케이션인 [Cowork](https://claude.com/product/cowork)용으로 설계되었지만 Claude Code에서도 동작합니다. 채용, 온보딩, 성과 관리, 정책 안내, 보상 분석을 돕습니다. 입력만으로 단독 동작하고 HRIS, ATS 및 다른 도구를 연결하면 더 강력해집니다.

## 설치

```bash
claude plugins add knowledge-work-plugins/human-resources
```

## 명령

Slash command로 호출하는 명시적 워크플로입니다.

| 명령 | 설명 |
|---|---|
| `/draft-offer` | 보상 세부사항, 시작일, 조건이 포함된 오퍼 레터 초안을 작성합니다. |
| `/onboarding` | 신규 입사자를 위한 온보딩 체크리스트와 첫 주 계획을 생성합니다. |
| `/performance-review` | 자기평가 프롬프트, 매니저 템플릿, 캘리브레이션 준비를 포함해 성과 평가를 구조화합니다. |
| `/policy-lookup` | PTO, 복리후생, 비용, 출장, 원격 근무 등 회사 정책을 찾고 설명합니다. |
| `/comp-analysis` | 벤치마킹, 밴드 배치, 지분 보상 갱신 모델링으로 보상 데이터를 분석합니다. |
| `/people-report` | 인원, 이탈, 다양성, 조직 건강도 보고서를 생성합니다. |

모든 명령은 context와 detail만 제공해도 **단독**으로 동작하며, MCP connector를 연결하면 더 강력해집니다.

## 스킬

관련 상황에서 Claude가 자동으로 사용하는 도메인 지식입니다.

| 스킬 | 설명 |
|---|---|
| `recruiting-pipeline` | Recruiting pipeline stage를 track하고 manage합니다. "recruiting update", "candidate pipeline", "how many candidates", "hiring status" 또는 sourcing, screening, interviewing, offer extension 논의 시 트리거됩니다. |
| `employee-handbook` | Company policy, benefit, procedure에 대한 질문에 답합니다. |
| `compensation-benchmarking` | Market data 기준으로 base, equity, total comp를 benchmark합니다. |
| `org-planning` | Headcount planning, org design, team structure optimization을 수행합니다. "org planning", "headcount plan", "team structure", "reorg", "who should we hire next" 또는 team size, reporting structure, organizational design을 고민할 때 트리거됩니다. |
| `people-analytics` | Workforce data를 분석해 attrition trend, engagement signal, diversity metric을 파악합니다. |
| `interview-prep` | 역량 기반 질문과 스코어카드가 포함된 구조화된 인터뷰 계획을 만듭니다. "interview plan for", "interview questions for", "how should we interview", "scorecard for" 또는 후보자 인터뷰 준비 요청에서 트리거됩니다. |

## 예시 워크플로

### 오퍼 초안 작성

```text
/draft-offer
```

Role, level, location, compensation detail을 알려주면 term, equity breakdown, benefit summary가 포함된 complete offer letter draft를 받습니다.

### 신규 입사자 온보딩

```text
/onboarding
```

New hire의 name, role, team, start date를 제공하면 comprehensive onboarding checklist, first-week calendar, tool access list, buddy assignment template을 받습니다.

### 성과 평가 준비

```text
/performance-review
```

Self-assessment, manager review, calibration용 template을 받습니다. 구체적이고 actionable하며 fair한 feedback 구조화를 돕습니다.

### 복리후생 이해

자연어로 물어보면 됩니다.

```text
우리 parental leave policy가 뭐야?
```

`employee-handbook` skill이 자동으로 trigger되어 connected knowledge base에서 답을 찾습니다.

### 보상 벤치마킹

```text
/comp-analysis
```

Comp data를 업로드하거나 band를 설명하면 market comparison, band placement analysis, adjustment recommendation을 받습니다.

## 단독 사용과 통합 사용

모든 명령과 스킬은 integration 없이도 동작합니다.

| 할 수 있는 일 | 단독 사용 | 통합 사용 시 |
|-----------------|------------|-------------------|
| Offer 초안 | Detail을 수동 제공 | HRIS, ATS로 auto-fill |
| Onboarding checklist | Process 설명 | HRIS, Knowledge base template |
| Performance review | 수동 입력 | HRIS review history |
| Policy lookup | Handbook content 붙여넣기 | Knowledge base |
| Comp analysis | CSV upload, band 설명 | Compensation data MCP |
| People report | Data upload | HRIS live data |

## MCP 통합

> 익숙하지 않은 placeholder가 보이거나 연결된 도구를 확인해야 한다면 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

더 풍부한 경험을 위해 tool을 연결하세요.

| 범주 | 예시 | 가능해지는 일 |
|---|---|---|
| **HRIS** | Workday, BambooHR, Rippling | Employee data, org structure, PTO balance |
| **ATS** | Greenhouse, Lever, Ashby | Candidate pipeline, interview schedule, offer tracking |
| **Compensation** | Pave, Radford | Market benchmark, comp band data |
| **Chat** | Slack, Teams | Team announcement, candidate coordination |
| **Calendar** | Google Calendar, Microsoft 365 | Interview scheduling, onboarding calendar |
| **Email** | Gmail, Microsoft 365 | Offer letter, candidate communication |

지원되는 integration 전체 목록은 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

## 설정

개인화를 위해 `settings.local.json` file을 만드세요.

- **Cowork**: Folder picker로 Cowork와 공유한 아무 folder에 저장합니다. Plugin이 자동으로 찾습니다.
- **Claude Code**: `human-resources/.claude/settings.local.json`에 저장합니다.

```json
{
  "company": "Your Company",
  "headquarters": "City, State",
  "employeeCount": 500,
  "benefits": {
    "healthInsurance": "Provider Name",
    "pto": "Unlimited / X days",
    "parentalLeave": "X weeks"
  },
  "compensation": {
    "currency": "USD",
    "equityType": "RSU / Options",
    "vestingSchedule": "4 years, 1 year cliff"
  }
}
```

설정되어 있지 않으면 플러그인이 이 정보를 대화형으로 물어봅니다.
