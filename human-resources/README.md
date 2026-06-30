# HR 플러그인

인사 운영 플러그인입니다. 주로 Anthropic의 에이전트형 데스크톱 애플리케이션인 [Cowork](https://claude.com/product/cowork)용으로 설계되었지만 Claude Code에서도 동작합니다. 채용, 온보딩, 성과 관리, 정책 안내, 보상 분석을 돕습니다. 입력만으로 단독 동작하고 HRIS, ATS 및 다른 도구를 연결하면 더 강력해집니다.

## 설치

```bash
claude plugins add knowledge-work-plugins/human-resources
```

## 명령

슬래시 명령으로 호출하는 명시적 워크플로입니다.

| 명령 | 설명 |
|---|---|
| `/draft-offer` | 보상 세부사항, 시작일, 조건이 포함된 오퍼 레터 초안을 작성합니다. |
| `/onboarding` | 신규 입사자를 위한 온보딩 체크리스트와 첫 주 계획을 생성합니다. |
| `/performance-review` | 자기평가 프롬프트, 매니저 템플릿, 캘리브레이션 준비를 포함해 성과 평가를 구조화합니다. |
| `/policy-lookup` | PTO, 복리후생, 비용, 출장, 원격 근무 등 회사 정책을 찾고 설명합니다. |
| `/comp-analysis` | 벤치마킹, 밴드 배치, 지분 보상 갱신 모델링으로 보상 데이터를 분석합니다. |
| `/people-report` | 인원, 이탈, 다양성, 조직 건강도 보고서를 생성합니다. |

모든 명령은 맥락과 세부사항만 제공해도 **단독**으로 동작하며, MCP 커넥터를 연결하면 더 강력해집니다.

## 스킬

관련 상황에서 Claude가 자동으로 사용하는 도메인 지식입니다.

| 스킬 | 설명 |
|---|---|
| `recruiting-pipeline` | 채용 파이프라인 단계를 추적하고 관리합니다. "recruiting update", "candidate pipeline", "how many candidates", "hiring status" 또는 소싱, 스크리닝, 면접, 오퍼 제안 논의 시 트리거됩니다. |
| `employee-handbook` | 회사 정책, 복리후생, 절차에 대한 질문에 답합니다. |
| `compensation-benchmarking` | 시장 데이터 기준으로 기본급, 지분 보상, 총보상을 벤치마킹합니다. |
| `org-planning` | 인력 계획, 조직 설계, 팀 구조 최적화를 수행합니다. "org planning", "headcount plan", "team structure", "reorg", "who should we hire next" 또는 팀 규모, 보고 체계, 조직 설계를 고민할 때 트리거됩니다. |
| `people-analytics` | 인력 데이터를 분석해 이탈 추세, 참여 신호, 다양성 지표를 파악합니다. |
| `interview-prep` | 역량 기반 질문과 스코어카드가 포함된 구조화된 인터뷰 계획을 만듭니다. "interview plan for", "interview questions for", "how should we interview", "scorecard for" 또는 후보자 인터뷰 준비 요청에서 트리거됩니다. |

## 예시 워크플로

### 오퍼 초안 작성

```text
/draft-offer
```

직무, 레벨, 근무지, 보상 세부사항을 알려주면 조건, 지분 보상 내역, 복리후생 요약이 포함된 완성형 오퍼 레터 초안을 받습니다.

### 신규 입사자 온보딩

```text
/onboarding
```

신규 입사자의 이름, 직무, 팀, 입사일을 제공하면 종합 온보딩 체크리스트, 첫 주 일정, 도구 접근 목록, 버디 배정 템플릿을 받습니다.

### 성과 평가 준비

```text
/performance-review
```

자기평가, 매니저 리뷰, 캘리브레이션용 템플릿을 받습니다. 구체적이고 실행 가능하며 공정한 피드백 구조화를 돕습니다.

### 복리후생 이해

자연어로 물어보면 됩니다.

```text
우리 육아휴직 정책이 뭐야?
```

`employee-handbook` 스킬이 자동으로 트리거되어 연결된 지식 베이스에서 답을 찾습니다.

### 보상 벤치마킹

```text
/comp-analysis
```

보상 데이터를 업로드하거나 밴드를 설명하면 시장 비교, 밴드 배치 분석, 조정 추천을 받습니다.

## 단독 사용과 통합 사용

모든 명령과 스킬은 통합 없이도 동작합니다.

| 할 수 있는 일 | 단독 사용 | 통합 사용 시 |
|-----------------|------------|-------------------|
| 오퍼 초안 | 세부사항을 수동 제공 | HRIS, ATS로 자동 입력 |
| 온보딩 체크리스트 | 프로세스 설명 | HRIS, 지식 베이스 템플릿 |
| 성과 평가 | 수동 입력 | HRIS 평가 이력 |
| 정책 조회 | 핸드북 내용 붙여넣기 | 지식 베이스 |
| 보상 분석 | CSV 업로드, 밴드 설명 | 보상 데이터 MCP |
| 인력 보고서 | 데이터 업로드 | HRIS 실시간 데이터 |

## MCP 통합

> 익숙하지 않은 플레이스홀더가 보이거나 연결된 도구를 확인해야 한다면 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

더 풍부한 경험을 위해 도구를 연결하세요.

| 범주 | 예시 | 가능해지는 일 |
|---|---|---|
| **HRIS** | Workday, BambooHR, Rippling | 직원 데이터, 조직 구조, PTO 잔여일 |
| **ATS** | Greenhouse, Lever, Ashby | 후보자 파이프라인, 면접 일정, 오퍼 추적 |
| **Compensation** | Pave, Radford | 시장 벤치마크, 보상 밴드 데이터 |
| **Chat** | Slack, Teams | 팀 공지, 후보자 커뮤니케이션 조율 |
| **Calendar** | Google Calendar, Microsoft 365 | 면접 일정 조율, 온보딩 캘린더 |
| **Email** | Gmail, Microsoft 365 | 오퍼 레터, 후보자 커뮤니케이션 |

지원되는 통합 전체 목록은 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

## 설정

개인화를 위해 `settings.local.json` 파일을 만드세요.

- **Cowork**: 폴더 선택기로 Cowork와 공유한 아무 폴더에 저장합니다. 플러그인이 자동으로 찾습니다.
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
