# Sales 플러그인

Anthropic의 에이전트 데스크톱 애플리케이션인 [Cowork](https://claude.com/product/cowork)을 위해 설계된 영업 생산성 플러그인이며, Claude Code에서도 사용할 수 있습니다. 잠재 고객 발굴, 아웃리치, 파이프라인 관리, 통화 준비, 거래 전략을 지원합니다. 어떤 영업 팀이든 사용할 수 있으며, 웹 검색과 사용자 입력만으로도 단독 사용이 가능하고, CRM, 이메일 등의 도구를 연결하면 더욱 강력해집니다.

## 설치

```bash
claude plugins add knowledge-work-plugins/sales
```

## 커맨드

슬래시 커맨드로 호출하는 명시적 워크플로우입니다:

| 커맨드 | 설명 |
|---|---|
| `/call-summary` | 통화 메모 또는 녹취록 처리 — 액션 아이템 추출, 후속 조치 이메일 작성, 내부 요약 생성 |
| `/forecast` | 가중 영업 예측 생성 — CSV 업로드 또는 파이프라인 설명, 할당량 설정, 전망 확인 |
| `/pipeline-review` | 파이프라인 상태 분석 — 거래 우선순위 지정, 리스크 플래그, 주간 실행 계획 수립 |

모든 커맨드는 **단독**으로 작동하며(메모 붙여넣기, CSV 업로드, 상황 설명), MCP 커넥터를 연결하면 **더욱 강력해집니다**.

## 스킬

관련 상황에서 Claude가 자동으로 활용하는 도메인 지식입니다:

| 스킬 | 설명 |
|---|---|
| `account-research` | 기업 또는 인물 리서치 — 기업 정보, 주요 담당자, 최신 뉴스, 채용 신호를 웹 검색으로 조사 |
| `call-prep` | 영업 통화 준비 — 어카운트 맥락, 참석자 조사, 제안 안건, 디스커버리 질문 |
| `daily-briefing` | 우선순위별 일일 영업 브리핑 — 미팅, 파이프라인 알림, 이메일 우선순위, 권장 조치 |
| `draft-outreach` | 리서치 기반 아웃리치 — 잠재 고객을 조사한 후 맞춤형 이메일 및 LinkedIn 메시지 작성 |
| `competitive-intelligence` | 경쟁사 조사 — 제품 비교, 가격 정보, 최근 출시, 차별화 매트릭스, 세일즈 토크 트랙 |
| `create-an-asset` | 맞춤형 영업 자산 생성 — 잠재 고객에 맞춘 랜딩 페이지, 덱, 원페이저, 워크플로우 데모 |

## 예시 워크플로우

### 통화 후

```
/call-summary
```

메모 또는 녹취록을 붙여넣으면 구조화된 요약, 담당자별 액션 아이템, 후속 조치 이메일 초안을 받을 수 있습니다. CRM이 연결되어 있으면 활동 기록 및 작업 생성을 제안합니다.

### 주간 예측

```
/forecast
```

CRM에서 내보낸 CSV를 업로드하거나 거래 내역을 붙여넣으세요. 할당량과 기간을 알려주시면 최선/예상/최악 시나리오에 따른 가중 예측, 확약 vs. 상향 분석, 갭 분석을 제공합니다.

### 파이프라인 리뷰

```
/pipeline-review
```

CSV를 업로드하거나 파이프라인을 설명하세요. 건강 점수, 거래 우선순위, 리스크 플래그(정체 거래, 마감일 경과, 단일 스레드), 주간 실행 계획을 받을 수 있습니다.

### 잠재 고객 리서치

자연스럽게 요청하세요:
```
내일 미팅 전에 Acme Corp에 대해 조사해줘
```

`account-research` 스킬이 자동으로 실행되어 기업 개요, 주요 담당자, 최신 뉴스, 권장 접근 방식을 제공합니다.

### 아웃리치 작성

```
TechStart의 Engineering VP에게 이메일을 작성해줘
```

`draft-outreach` 스킬이 먼저 잠재 고객을 조사한 후 다양한 접근 각도의 맞춤형 아웃리치를 생성합니다.

### 경쟁 정보

```
TechStart와 Competitor X를 비교해줘
```

`competitive-intelligence` 스킬이 양사를 조사하고 토크 트랙이 포함된 차별화 매트릭스를 구축합니다.

## 단독 사용 + 강화 모드

모든 커맨드와 스킬은 통합 도구 없이도 작동합니다:

| 가능한 작업 | 단독 사용 | 강화 연동 |
|-------------|-----------|-----------|
| 통화 메모 처리 | 메모/녹취록 붙여넣기 | Transcripts MCP (예: Gong, Fireflies) |
| 파이프라인 예측 | CSV 업로드, 거래 붙여넣기 | CRM MCP |
| 파이프라인 리뷰 | CSV 업로드, 거래 설명 | CRM MCP |
| 잠재 고객 조사 | 웹 검색 | Enrichment MCP (예: Clay, ZoomInfo) |
| 통화 준비 | 미팅 설명 | CRM, Email, Calendar MCP |
| 아웃리치 작성 | 웹 검색 + 사용자 맥락 | CRM, Email MCP |
| 경쟁 정보 | 웹 검색 | CRM (성패 데이터), Docs (배틀카드) |

## MCP 통합

> 익숙하지 않은 플레이스홀더가 보이거나 연결된 도구를 확인하려면 [CONNECTORS.md](CONNECTORS.md)를 참조하세요.

도구를 연결하면 더 풍부한 경험을 제공합니다:

| 카테고리 | 예시 | 활성화되는 기능 |
|---|---|---|
| **CRM** | HubSpot, Close | 파이프라인 데이터, 어카운트 이력, 연락처 기록 |
| **녹취록** | Fireflies, Gong, Chorus | 통화 녹음, 녹취록, 핵심 순간 |
| **데이터 보강** | Clay, ZoomInfo, Apollo | 기업 및 연락처 데이터 보강 |
| **채팅** | Slack, Teams | 내부 논의, 동료 인텔리전스 |

이메일, 캘린더 및 추가 CRM 옵션을 포함한 전체 지원 통합 목록은 [CONNECTORS.md](CONNECTORS.md)를 참조하세요.

## 설정

개인 설정을 위해 `sales/.claude/settings.local.json`에 로컬 설정 파일을 생성하세요:

```json
{
  "name": "Your Name",
  "title": "Account Executive",
  "company": "Your Company",
  "quota": {
    "annual": 1000000,
    "quarterly": 250000
  },
  "product": {
    "name": "Your Product",
    "value_props": [
      "Key value proposition 1",
      "Key value proposition 2"
    ],
    "competitors": [
      "Competitor A",
      "Competitor B"
    ]
  }
}
```

설정이 구성되지 않은 경우 플러그인이 대화형으로 해당 정보를 요청합니다.
