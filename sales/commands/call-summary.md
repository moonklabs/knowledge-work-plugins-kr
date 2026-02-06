---
description: 콜 노트 또는 녹취를 처리해 액션 아이템 추출, 후속 이메일 초안, 내부 요약 생성
argument-hint: "<call notes or transcript>"
---

# /call-summary

> 낯선 플레이스홀더가 보이거나 어떤 도구가 연결되어 있는지 확인하려면 [CONNECTORS.md](../CONNECTORS.md)를 보세요.

콜 노트 또는 녹취를 처리해 액션 아이템을 추출하고, 후속 커뮤니케이션을 초안 작성하며, 기록을 업데이트합니다.

## 사용법

```
/call-summary
```

그 다음 노트나 녹취를 붙여넣거나 녹음을 업로드하세요.

---

## 동작 방식

```
┌─────────────────────────────────────────────────────────────────┐
│                      CALL SUMMARY                                │
├─────────────────────────────────────────────────────────────────┤
│  STANDALONE (always works)                                       │
│  ✓ 콜 노트 또는 녹취를 붙여넣기                                 │
│  ✓ 핵심 논의 포인트 및 결정 사항 추출                           │
│  ✓ 담당자와 마감일이 포함된 액션 아이템 식별                     │
│  ✓ 이의사항, 우려, 열린 질문 노출                               │
│  ✓ 고객용 후속 이메일 초안 작성                                 │
│  ✓ 팀을 위한 내부 요약 생성                                     │
├─────────────────────────────────────────────────────────────────┤
│  SUPERCHARGED (when you connect your tools)                      │
│  + Transcripts: 녹음을 자동으로 가져오기 (예: Gong, Fireflies)    │
│  + CRM: 기회 업데이트, 활동 로그 기록, 작업 생성                │
│  + Email: 초안에서 바로 후속 메일 발송                           │
│  + Calendar: 미팅 연결, 참석자 컨텍스트 가져오기                 │
└─────────────────────────────────────────────────────────────────┘
```

---

## 필요한 정보

**옵션 1: 노트 붙여넣기**
어떤 형태든 붙여넣어 주세요. 불릿, 러프 노트, 생각나는 대로 적은 글 모두 가능합니다. 제가 구조화합니다.

**옵션 2: 녹취 붙여넣기**
화상 회의 도구(예: Zoom, Teams) 또는 대화 인텔 도구(예: Gong, Fireflies)에서 전체 녹취가 있다면 붙여넣어 주세요. 핵심 구간을 추출합니다.

**옵션 3: 콜 설명**
무슨 일이 있었는지 간단히 설명해 주세요: "Acme Corp와 디스커버리 콜을 했습니다. VP Eng와 CTO를 만났고, Competitor X와 비교 중입니다. 주요 우려는 통합 일정입니다."

---

## 출력

### Internal Summary
```markdown
## Call Summary: [Company] — [Date]

**Attendees:** [Names and titles]
**Call Type:** [Discovery / Demo / Negotiation / Check-in]
**Duration:** [If known]

### Key Discussion Points
1. [Topic] — [What was discussed, decisions made]
2. [Topic] — [Summary]

### Customer Priorities
- [Priority 1 they expressed]
- [Priority 2]

### Objections / Concerns Raised
- [Concern] — [How you addressed it / status]

### Competitive Intel
- [Any competitor mentions, what was said]

### Action Items
| Owner | Action | Due |
|-------|--------|-----|
| [You] | [Task] | [Date] |
| [Customer] | [Task] | [Date] |

### Next Steps
- [Agreed next step with timeline]

### Deal Impact
- [How this call affects the opportunity — stage change, risk, acceleration]
```

### Customer Follow-Up Email
```
Subject: [Meeting recap + next steps]

Hi [Name],

Thank you for taking the time to meet today...

[Key points discussed]

[Commitments you made]

[Clear next step with timeline]

Best,
[You]
```

---

## 이메일 스타일 가이드라인

고객용 이메일을 작성할 때:

1. **간결하지만 정보는 충분히** — 바로 요점을 전달하세요. 고객은 바쁩니다.
2. **마크다운 금지** — 별표, 굵게 등 마크다운 문법을 쓰지 마세요. 어떤 이메일 클라이언트에서도 자연스럽게 보이는 일반 텍스트로 작성합니다.
3. **단순한 구조** — 짧은 문단, 섹션 사이 줄바꿈. 고객의 이메일 클라이언트가 렌더링하지 않는다면 헤더/불릿을 사용하지 않습니다.
4. **읽기 쉽게** — 항목 나열 시 화려한 포맷 대신 일반 대시나 숫자를 사용합니다.

**좋은 예:**
```
Here's what we discussed:
- Quote for 20 seats at $480/seat/year
- W9 and supplier onboarding docs
- Point of contact for the contract
```

**나쁜 예:**
```
**What You Need from Us:**
- Quote for 20 seats at $480/seat/year
```

---

## 커넥터 연결 시

**Transcripts 연결됨(예: Gong, Fireflies):**
- 콜을 자동으로 찾습니다
- 전체 녹취를 가져옵니다
- 플랫폼이 표시한 핵심 구간을 추출합니다

**CRM 연결됨:**
- 기회 단계 업데이트를 제안합니다
- 콜을 활동 로그로 기록합니다
- 액션 아이템을 작업으로 생성합니다
- 다음 단계 필드를 업데이트합니다

**Email 연결됨:**
- ~~email에서 초안을 만들도록 제안합니다
- 승인 시 바로 발송할 수 있습니다

---

## 팁

1. **상세할수록 좋습니다** — 러프 노트도 도움이 됩니다. "X가 걱정되는 것 같았다" 같은 컨텍스트가 유용합니다.
2. **참석자를 알려주세요** — 요약 구조화와 액션 아이템 할당에 도움이 됩니다.
3. **중요한 포인트 표시** — 중요했던 내용은 "가장 큰 이슈는..."처럼 알려 주세요.
4. **딜 스테이지를 알려주세요** — 후속 톤과 다음 단계를 맞추는 데 도움이 됩니다.
