# Operations 플러그인

비즈니스 운영 플러그인입니다. 주로 Anthropic의 agentic 데스크톱 애플리케이션인 [Cowork](https://claude.com/product/cowork)용으로 설계되었지만 Claude Code에서도 동작합니다. 벤더 관리, 프로세스 문서화, 변경 관리, 수용량 계획, 컴플라이언스 추적, 리소스 계획을 돕습니다. 입력만으로 단독 동작하고 ITSM, 프로젝트 추적기 및 다른 도구를 연결하면 더 강력해집니다.

## 설치

```bash
claude plugins add knowledge-work-plugins/operations
```

## 명령

슬래시 명령으로 호출하는 명시적 워크플로입니다.

| 명령 | 설명 |
|---|---|
| `/vendor-review` | 비용 분석, 위험 평가, 계약 요약, 갱신 추천으로 벤더를 평가합니다. |
| `/process-doc` | 흐름도, RACI 매트릭스, SOP, 런북으로 비즈니스 프로세스를 문서화합니다. |
| `/change-request` | 영향 분석, 롤백 계획, 승인 라우팅이 포함된 변경 관리 요청을 작성합니다. |
| `/capacity-plan` | 업무량 분석, 인원 모델링, 활용률 예측으로 리소스 수용량을 계획합니다. |
| `/status-report` | 리더십용 프로젝트 업데이트, KPI, 위험, 액션 아이템이 포함된 상태 보고서를 생성합니다. |
| `/runbook` | 반복 작업을 위한 단계별 운영 런북을 생성하거나 업데이트합니다. |

모든 명령은 맥락과 세부 정보만 제공해도 **단독**으로 동작하며, MCP connector를 연결하면 더 강력해집니다.

## 스킬

관련 상황에서 Claude가 자동으로 사용하는 도메인 지식입니다.

| 스킬 | 설명 |
|---|---|
| `vendor-management` | Contract, performance, risk 기준으로 vendor relationship을 evaluate, compare, manage합니다. |
| `process-optimization` | 비즈니스 프로세스를 분석하고 개선합니다. "this process is slow", "how can we improve", "streamline this workflow", "too many steps", "bottleneck" 또는 사용자가 고치고 싶은 비효율적인 프로세스를 설명할 때 트리거됩니다. |
| `change-management` | Communication, training, adoption을 포함해 organizational 또는 technical change를 계획하고 실행합니다. |
| `risk-assessment` | 운영 위험을 식별, 평가, 완화합니다. "what are the risks", "risk assessment", "risk register", "what could go wrong" 또는 프로젝트, 벤더, 프로세스, 의사결정 관련 위험을 평가할 때 트리거됩니다. |
| `compliance-tracking` | Compliance requirement와 audit readiness를 track합니다. "compliance", "audit prep", "SOC 2", "ISO 27001", "GDPR", "regulatory requirement" 또는 compliance activity tracking/preparation/documentation에 도움이 필요할 때 트리거됩니다. |
| `resource-planning` | 수용량, 활용률, 예측, 예산을 기준으로 resource allocation을 계획하고 optimize합니다. |

## 예시 워크플로

### 벤더 평가

```text
/vendor-review
```

벤더 이름과 계약 세부 정보를 제공하거나 제안서를 업로드하면 비용 분석, 위험 표시, 추천 조치가 포함된 구조화된 평가를 받습니다.

### 프로세스 문서화

```text
/process-doc employee offboarding
```

프로세스를 설명하거나 단계별로 알려주면 flowchart, RACI matrix, 단계별 절차가 포함된 완성된 SOP를 받습니다.

### 변경 요청 제출

```text
/change-request
```

변경 사항을 설명하면 approval에 바로 쓸 수 있는 impact analysis, risk assessment, rollback plan, communication template을 받습니다.

### 수용량 계획

```text
/capacity-plan
```

팀 데이터를 업로드하거나 resource를 설명하면 활용률 분석, 병목 식별, headcount recommendation을 받습니다.

### 리더십 상태 보고서

```text
/status-report
```

연결된 도구에서 update를 가져오거나 필요한 input을 물은 뒤 KPI, risk, next step이 포함된 다듬어진 상태 보고서를 생성합니다.

### 런북 작성

```text
/runbook monthly close process
```

프로세스를 한 번 설명하면 checklist, troubleshooting, escalation path가 포함된 반복 가능한 runbook으로 문서화합니다.

## 단독 사용과 통합 사용

모든 명령과 스킬은 통합 도구 없이도 동작합니다.

| 할 수 있는 일 | 단독 사용 | 통합 사용 시 |
|-----------------|------------|-------------------|
| Vendor review | Detail 제공, proposal upload | Procurement, Knowledge base |
| Process documentation | Process 설명 | Knowledge base(existing docs) |
| Change request | Change 설명 | ITSM, Project tracker |
| Capacity planning | Data upload, team 설명 | Project tracker(workload data) |
| Status report | Update 수동 제공 | Project tracker, Chat, Calendar |
| Runbook | Process walkthrough | Knowledge base, ITSM |

## MCP 통합

> 익숙하지 않은 placeholder가 보이거나 연결된 도구를 확인해야 한다면 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

더 풍부한 경험을 위해 도구를 연결하세요.

| 범주 | 예시 | 가능해지는 일 |
|---|---|---|
| **ITSM** | ServiceNow, Zendesk | Ticket management, change request, incident tracking |
| **Project tracker** | Asana, Jira, monday.com | Project status, resource allocation, task tracking |
| **Knowledge base** | Notion, Confluence | Process doc, runbook, policy |
| **Chat** | Slack, Teams | Team coordination, approval, status update |
| **Calendar** | Google Calendar, Microsoft 365 | Meeting scheduling, deadline tracking |
| **Email** | Gmail, Microsoft 365 | Vendor communication, approval |

지원되는 통합 전체 목록은 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

## 설정

개인화를 위해 `operations/.claude/settings.local.json`에 로컬 설정 파일을 만드세요.

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
