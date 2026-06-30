# 소규모 사업 플러그인

[Cowork](https://claude.com/product/cowork)용으로 미리 구성된 소규모 사업 워크플로 모음입니다. Cowork는 Anthropic의 에이전트형 데스크톱 애플리케이션이며, 이 플러그인은 Claude Code에서도 동작합니다. 설치하면 15개의 구성요소 스킬, 15개의 바로 쓸 수 있는 워크플로, 자연어 요청을 이해하는 라우터를 사용할 수 있습니다.

무엇을 외울 필요가 없습니다. Claude에게 필요한 일을 그대로 말하면 됩니다. 예: "급여를 지급할 수 있을지 걱정돼", "고객이 화가 났어", "가격을 얼마로 잡아야 할까?". Claude가 적절한 워크플로를 찾고 단계별로 안내합니다. 모든 워크플로는 작업 전에 멈추므로, 사업주의 승인 없이 실행되는 일은 없습니다.

> **중요**: 이 플러그인은 소규모 사업 워크플로를 보조하지만 재무, 세무, 법무, HR 자문을 제공하지 않습니다. 모든 출력 결과는 사용 전 직접 검토해야 하며, 필요한 경우 자격 있는 전문가의 검토를 받아야 합니다.

## 설치

### Cowork

[claude.com/plugins](https://claude.com/plugins/)에서 설치합니다.

### Claude Code

```bash
claude plugin marketplace add anthropics/knowledge-work-plugins
claude plugin install small-business@knowledge-work-plugins
```

설치가 끝나면 **"set me up"**이라고 말해 `smb-onboard` 스킬을 실행하세요. Claude가 사업, 고민, 이미 사용하는 도구를 이해하도록 도와줍니다.

## 연결하면 좋은 도구

`/smb-onboard`를 실행하거나 Claude에게 "set me up"이라고 말하세요.

**핵심 도구** (가장 좋은 경험을 위해 먼저 연결):
- **QuickBooks** — 현금 예측, 마진, 월말 마감, 세금 준비 등 재무 워크플로를 구동합니다
- **PayPal** — 거래 데이터, 인보이스, 분쟁, 환불
- **HubSpot** — CRM, 리드, 캠페인, 고객 지원 티켓

**마케팅 및 커뮤니케이션:**
- **Canva** — 브랜드에 맞는 소셜 및 이메일 자산 생성
- **Gmail / Outlook** — 이메일 초안, 티켓 처리, 계약 검토
- **Google Calendar / Outlook Calendar** — 미팅 준비, 통화 시간 확보, 주간 약속 확인
- **Slack** — 브리프 전달 및 알림

**선택 도구** (연결하면 더 깊은 분석 가능):
- **Stripe** — 결제 및 구독 데이터
- **Square** — POS 거래 데이터
- **Google Drive / OneDrive** — 파일 저장소와 템플릿
- **DocuSign** — 대기 중인 envelope의 계약 검토

처음부터 전부 연결할 필요는 없습니다. 한두 개만 연결해도 바로 가치가 보이며, 다른 도구를 연결하면 더 많은 기능이 열릴 때 플러그인이 알려줍니다.

## 작동 방식

세 계층이 함께 동작합니다.

1. **스킬** — 구성요소입니다. 각 스킬은 현금 예측, 리드 점수화, 인보이스 알림 초안처럼 한 가지 일을 잘 처리합니다. 총 15개가 있습니다.

2. **명령** — 워크플로입니다. 명령은 여러 스킬을 다단계 절차로 연결하고, 실행 전 승인 체크포인트를 둡니다. 총 15개가 있습니다.

3. **라우터** — 첫 진입점입니다. Claude에게 자연어로 말하면 라우터가 적합한 워크플로를 찾아 안내합니다. 명령 이름을 외울 필요가 없습니다.

## 전체 15개 명령

명령은 여러 스킬을 연결하는 워크플로입니다. 각 명령은 작업 전에 승인 체크포인트에서 멈춥니다.

### 자금 및 재무

| 명령 | 수행 내용 | 이렇게 말해보세요 | 사용 스킬 | 필수 | 선택 |
|---|---|---|---|---|---|
| `/plan-payroll` | 급여 지급 가능 여부를 판단할 수 있도록 현금을 예측하고 연체 인보이스 회수를 준비합니다. | "급여 줄 수 있을까", "현금이 빠듯해", "누가 아직 돈을 안 냈지" | cash-flow-snapshot, invoice-chase | QuickBooks | PayPal, Stripe, Square, Mail |
| `/month-heads-up` | 30일 현금 전망과 조기 위험 신호를 제공합니다. | "다음 달 현금흐름 어때", "현금 예측해줘", "버틸 수 있는 기간이 궁금해" | cash-flow-snapshot | QuickBooks | PayPal |
| `/close-month` | 월말 마감: 대사, 누락 표시, P&L 작성, 마감 패킷 내보내기를 수행합니다. | "장부 마감해줘", "월말 마감", "대사해줘" | month-end-prep | QuickBooks | PayPal, Stripe, Square |
| `/price-check` | 제품별 마진 표와 세 가지 가격 시나리오를 제공합니다. | "마진이 어떻게 돼", "가격을 올려야 할까", "단위당 비용이 궁금해" | margin-analyzer | QuickBooks | PayPal |
| `/tax-prep` | 회계사에게 넘길 세금 준비 자료를 만듭니다(분기 추정세 또는 연말 1099). | "세금 준비", "추정세", "1099", "회계사에게 줄 자료가 필요해" | tax-season-organizer | QuickBooks | PayPal, Stripe |

### 영업 및 마케팅

| 명령 | 수행 내용 | 이렇게 말해보세요 | 사용 스킬 | 필수 | 선택 |
|---|---|---|---|---|---|
| `/call-list` | 오늘 전화할 상위 5개 리드와 대화 포인트, 캘린더 시간을 준비합니다. | "누구에게 전화하지", "뜨거운 리드 있어", "파이프라인 확인해줘" | lead-triage | HubSpot | Mail, Google Calendar |
| `/run-campaign` | 영업 분석 -> 콘텐츠 브리프 -> Canva 자산 -> HubSpot 발송까지 캠페인을 끝까지 준비합니다. | "캠페인 실행해줘", "매출이 떨어졌어", "고객이 더 필요해" | content-strategy, canva-creator, lead-triage | HubSpot, Canva | QuickBooks, PayPal |
| `/sales-brief` | 상위/하위 판매 상품과 2주 콘텐츠 브리프를 제공합니다. | "뭐가 잘 팔려", "무엇을 홍보하지" | content-strategy | QuickBooks 또는 PayPal | HubSpot |

### 고객 및 운영

| 명령 | 수행 내용 | 이렇게 말해보세요 | 사용 스킬 | 필수 | 선택 |
|---|---|---|---|---|---|
| `/customer-pulse-check` | 고객 피드백 주제와 답변 템플릿을 정리합니다. | "고객들이 뭐라고 해", "불만 사항", "리뷰 확인해줘" | customer-pulse, ticket-deflector | PayPal 또는 HubSpot | -- |
| `/handle-complaint` | 고객 불만 처리: 맥락 확인, 답변 초안, 운영 개선안 제안을 끝까지 수행합니다. | "고객이 화났어", "이 불만 처리해줘", "화난 이메일이 왔어" | ticket-deflector, customer-pulse | -- (붙여 넣은 텍스트로 동작) | Gmail, HubSpot, PayPal |
| `/crm-cleanup` | HubSpot 정리: 오래된 거래, 중복, 누락 필드를 찾고 승인한 항목만 수정합니다. | "CRM 정리해줘", "HubSpot이 엉망이야", "오래된 거래 확인" | crm-maintenance | HubSpot | -- |
| `/review-contract` | 쉬운 말로 계약을 검토하고 위험 신호와 심각도를 표시합니다. | "이 계약 검토해줘", "NDA", "서명해도 될까" | contract-review | -- (파일 업로드로 동작) | DocuSign |

### 비즈니스 인텔리전스

| 명령 | 수행 내용 | 이렇게 말해보세요 | 사용 스킬 | 필수 | 선택 |
|---|---|---|---|---|---|
| `/monday-brief` | 월요일 아침 브리핑: 현금, 매출, 파이프라인, 이번 주 일정, 상위 3개 할 일을 정리합니다. | "월요일 브리프", "오늘 뭐 해야 해", "이번 주 시작 브리핑" | business-pulse | -- (가능한 범위로 축소 동작) | QuickBooks, PayPal, HubSpot, Calendar, Gmail, Slack |
| `/friday-brief` | 금요일 주간 마감 pulse: 전주 대비 매출, 성과, 주의할 항목을 정리합니다. | "이번 주 마감", "이번 주 어땠어", "금요일 recap" | business-pulse | PayPal 또는 HubSpot | -- |
| `/quarterly-review` | 매출, 마진, 고객 건전성, 기회, 위험을 담은 전체 QBR 내러티브를 만듭니다. | "분기 리뷰", "이사회 덱", "QBR" | business-pulse | QuickBooks | PayPal, HubSpot |

## 전체 15개 스킬

스킬은 원자적 구성요소입니다. 각 스킬은 한 가지 일을 잘 처리합니다.

### 자금 및 재무

| 스킬 | 수행 내용 | 이렇게 말해보세요 | 필수 | 선택 |
|---|---|---|---|---|
| **cash-flow-snapshot** | 신뢰 구간과 이름 붙은 위험 신호가 포함된 30/60/90일 현금 예측을 만들고 채팅 요약과 XLSX를 제공합니다. | "현금흐름 예측해줘", "급여 줄 수 있을까", "런웨이", "현금 위기" | QuickBooks, PayPal, Stripe, Square 중 하나 | 나머지는 보조 소스 |
| **invoice-chase** | 각 고객의 결제 이력과 톤에 맞춘 연체 인보이스 알림 초안을 작성하고 승인 후 PayPal로 보냅니다. | "누가 돈을 안 냈지", "연체 인보이스", "미납 후속 연락" | QuickBooks | PayPal, Stripe, Gmail |
| **margin-analyzer** | 제품 또는 서비스별 단위 경제성을 인플레이션 기준과 비교하고 세 가지 가격 시나리오를 제공합니다. | "마진이 어떻게 돼", "가격을 올릴까", "비용이 이익을 깎고 있어", "얼마 받아야 해" | QuickBooks | PayPal, Square, CSV 업로드 |
| **month-end-prep** | 월말 마감: 결제 처리기관과 대사하고 누락을 표시하며 P&L 내러티브와 마감 패킷을 내보냅니다. | "이번 달 마감", "대사", "P&L", "매출이 왜 바뀌었지" | QuickBooks | PayPal, Stripe, Square |
| **tax-season-organizer** | 분기 추정세 계산 또는 연말 1099-NEC 준비와 회계사용 인수인계 패킷을 만듭니다. | "분기 세금", "추정세 납부", "1099", "1099-NEC", "연말 세금 준비" | QuickBooks | PayPal, Stripe |

### 영업 및 마케팅

| 스킬 | 수행 내용 | 이렇게 말해보세요 | 필수 | 선택 |
|---|---|---|---|---|
| **lead-triage** | 참여도, 적합도, 긴급도를 기준으로 HubSpot 리드를 점수화하고 대화 포인트가 포함된 전화 목록을 만듭니다. | "리드 우선순위", "누구부터 전화하지", "파이프라인" | HubSpot | Gmail, Google Calendar |
| **content-strategy** | 영업 데이터를 분석해 잘 팔리는 항목과 부진 항목을 찾고 우선순위가 있는 30일 콘텐츠 브리프를 만듭니다. | "무엇을 올리지", "콘텐츠 계획", "뭐가 팔려", "무엇을 홍보하지" | QuickBooks 또는 PayPal | Square |
| **canva-creator** | 콘텐츠 브리프를 받아 게시 캘린더, Canva 자산, 캡션 문구, HubSpot 스테이징까지 캠페인을 실행합니다. | "콘텐츠 만들어줘", "게시물 생성해줘", "자산 만들어줘", "캠페인으로 바꿔줘" | Canva, HubSpot | -- |

### 고객 및 운영

| 스킬 | 수행 내용 | 이렇게 말해보세요 | 필수 | 선택 |
|---|---|---|---|---|
| **customer-pulse** | 분쟁, 티켓, 이메일 감정, 리뷰를 모아 주제 보고서와 "이번 주 해야 할 세 가지" 목록을 만듭니다. | "고객 반응이 어때", "사람들이 뭐라고 해", "분쟁", "리뷰 분석" | -- (가능한 범위로 축소 동작) | PayPal, HubSpot, Gmail |
| **ticket-deflector** | 고객 이메일이나 티켓을 읽고 주문/환불 상태를 확인한 뒤 톤에 맞는 답변 초안을 작성합니다. 승인 시 PayPal 환불도 가능합니다. | "답장 초안 써줘", "이 고객에게 답해줘", "내 주문 어디 있냐고 해", "환불 원해" | PayPal, HubSpot, Mail | Intercom, Square |
| **crm-maintenance** | HubSpot을 최신 상태로 유지합니다. 연락처와 거래를 생성/수정하고 통화와 메모를 기록하며 오래된 레코드를 표시합니다. | "CRM 업데이트", "통화 기록", "HubSpot 정리", "거래에 맥락 추가" | HubSpot | Gmail, Google Calendar |
| **contract-review** | 쉬운 말로 계약을 검토하고 위험 신호, 심각도, 제안 수정안이 포함된 redline DOCX를 만듭니다. | "이 계약 검토해줘", "내가 뭘 서명하는 거야", "우려 사항 표시", "지급 조건 확인" | -- (파일 업로드로 동작) | Gmail, DocuSign |

### 채용

| 스킬 | 수행 내용 | 이렇게 말해보세요 | 필수 | 선택 |
|---|---|---|---|---|
| **job-post-builder** | 채용 공고, 점수표가 있는 구조화된 면접 가이드, 오퍼 레터 템플릿까지 전체 채용 패킷을 만듭니다. | "채용 도와줘", "채용 공고 써줘", "직무 설명", "오픈 포지션", "면접 질문", "오퍼 레터 초안" | -- (단독 동작) | DocuSign, Google Drive |

### 비즈니스 인텔리전스 및 온보딩

| 스킬 | 수행 내용 | 이렇게 말해보세요 | 필수 | 선택 |
|---|---|---|---|---|
| **business-pulse** | 현금, 매출, 파이프라인, 약속, 주의 목록, 오늘 가장 중요한 한 가지를 담은 1페이지 사업 스냅샷입니다. | "사업 상태 어때", "스냅샷", "주간 요약", "따라잡게 정리해줘" | -- (가능한 범위로 축소 동작) | QuickBooks, PayPal, HubSpot, Google Calendar, Gmail, Slack |
| **smb-onboard** | 도구 연결, 데모 절차 실행, 사업 맥락 저장, 주간 체크인 주기 설정을 안내합니다. | "set me up", "setup", "시작하자", "설정 도와줘", "처음이야", "뭘 할 수 있어" | -- | 모든 커넥터 |

## 커스터마이징

이 워크플로들은 일반적인 출발점입니다. 실제 사업 방식에 맞게 커스터마이즈하면 훨씬 유용해집니다.

- **사업 컨텍스트 추가** — 업종, 제품, 고객, 프로세스를 스킬 파일에 넣어 Claude가 맥락을 이해하게 합니다.
- **임계값 조정** — `business-pulse`와 `cash-flow-snapshot`의 알림 임계값을 규모에 맞게 조정합니다.
- **커넥터 교체** — 실제 사용하는 도구를 스킬이 바라보게 설정합니다.
