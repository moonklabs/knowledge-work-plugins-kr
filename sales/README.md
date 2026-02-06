# Sales 플러그인

[Cowork](https://claude.com/product/cowork)를 위한 세일즈 생산성 플러그인입니다. Anthropic의 에이전트형 데스크톱 앱인 Cowork에 맞춰 설계했지만 Claude Code에서도 동작합니다. 잠재 고객 발굴, 아웃리치, 파이프라인 관리, 콜 준비, 딜 전략 수립에 도움을 줍니다. 웹 검색과 사용자 입력만으로도 동작하며, CRM, 이메일 등 도구를 연결하면 성능이 크게 강화됩니다.

## 설치

```bash
claude plugins add knowledge-work-plugins/sales
```

## 커맨드

슬래시 커맨드로 실행하는 명시적 워크플로입니다.

| Command | 설명 |
|---|---|
| `/call-summary` | 콜 노트나 녹취를 처리해 액션 아이템 추출, 후속 메일 초안, 내부 요약 생성 |
| `/forecast` | 가중치 기반 세일즈 예측 생성: CSV 업로드 또는 파이프라인 설명, 쿼터 설정, 전망치 산출 |
| `/pipeline-review` | 파이프라인 건전성 분석: 딜 우선순위 지정, 리스크 플래그, 주간 액션 플랜 |

모든 커맨드는 **단독으로** 동작하며(노트 붙여넣기, CSV 업로드, 상황 설명), MCP 커넥터를 연결하면 **더 강력하게** 동작합니다.

## 스킬

관련 상황에서 자동으로 사용하는 도메인 지식입니다.

| Skill | 설명 |
|---|---|
| `account-research` | 회사 또는 인물 조사: 웹 검색으로 기업 정보, 주요 연락처, 최근 뉴스, 채용 신호 확인 |
| `call-prep` | 세일즈 콜 준비: 계정 컨텍스트, 참석자 조사, 추천 어젠다, 디스커버리 질문 |
| `daily-briefing` | 우선순위 기반 데일리 브리핑: 미팅, 파이프라인 알림, 이메일 우선순위, 추천 액션 |
| `draft-outreach` | 리서치 우선 아웃리치: 먼저 대상 조사 후 개인화 이메일/LinkedIn 메시지 작성 |
| `competitive-intelligence` | 경쟁사 조사: 제품 비교, 가격 정보, 최근 릴리스, 차별화 매트릭스, 세일즈 토크 트랙 |
| `create-an-asset` | 맞춤형 세일즈 에셋 생성: 랜딩 페이지, 덱, 원페이저, 워크플로 데모 |

## 예시 워크플로

### 콜 이후

```
/call-summary
```

노트나 녹취를 붙여넣으세요. 구조화된 요약, 담당자 포함 액션 아이템, 후속 이메일 초안을 제공합니다. CRM이 연결되어 있으면 활동 로그 및 작업 생성도 제안합니다.

### 주간 예측

```
/forecast
```

CRM에서 CSV를 내보내 업로드(또는 딜을 직접 붙여넣기)하세요. 쿼터와 기간을 알려주면 최선/가능/최악 시나리오 기반 가중 예측, 커밋 vs. 업사이드 분석, 갭 분석을 제공합니다.

### 파이프라인 리뷰

```
/pipeline-review
```

CSV를 업로드하거나 파이프라인을 설명하세요. 건강도 점수, 딜 우선순위, 리스크 플래그(오래된 딜, 지난 클로즈 날짜, 단일 접점 등), 주간 액션 플랜을 제공합니다.

### 잠재 고객 조사

자연스럽게 요청하세요.
```
내일 콜 전에 Acme Corp 조사해 주세요
```

`account-research` 스킬이 자동으로 실행되어 회사 개요, 주요 연락처, 최근 뉴스, 추천 접근 방법을 제공합니다.

### 아웃리치 초안

```
TechStart 엔지니어링 VP에게 보낼 이메일 초안 작성해 주세요
```

`draft-outreach` 스킬이 먼저 대상을 조사한 뒤, 여러 각도의 개인화 아웃리치를 생성합니다.

### 경쟁사 인텔

```
Competitor X와 비교하면 어떤가요?
```

`competitive-intelligence` 스킬이 두 회사의 정보를 조사해 차별화 매트릭스와 토크 트랙을 구성합니다.

## Standalone + Supercharged

모든 커맨드와 스킬은 통합 없이도 동작합니다.

| 가능한 작업 | Standalone | Supercharged With |
|-----------------|------------|-------------------|
| 콜 노트 처리 | 노트/녹취 붙여넣기 | Transcripts MCP (예: Gong, Fireflies) |
| 파이프라인 예측 | CSV 업로드, 딜 붙여넣기 | CRM MCP |
| 파이프라인 리뷰 | CSV 업로드, 딜 설명 | CRM MCP |
| 잠재 고객 조사 | 웹 검색 | Enrichment MCP (예: Clay, ZoomInfo) |
| 콜 준비 | 미팅 설명 | CRM, Email, Calendar MCPs |
| 아웃리치 작성 | 웹 검색 + 컨텍스트 | CRM, Email MCPs |
| 경쟁 인텔 | 웹 검색 | CRM(Win/Loss 데이터), Docs(배틀카드) |

## MCP 통합

> 낯선 플레이스홀더가 보이거나 어떤 도구가 연결되어 있는지 확인하려면 [CONNECTORS.md](CONNECTORS.md)를 보세요.

더 풍부한 경험을 위해 도구를 연결하세요.

| 카테고리 | 예시 | 지원 내용 |
|---|---|---|
| **CRM** | HubSpot, Close | 파이프라인 데이터, 계정 히스토리, 연락처 기록 |
| **Transcripts** | Fireflies, Gong, Chorus | 통화 녹음, 녹취록, 핵심 구간 |
| **Enrichment** | Clay, ZoomInfo, Apollo | 회사/연락처 데이터 보강 |
| **Chat** | Slack, Teams | 내부 논의 및 동료 인텔 |

이메일, 캘린더, 추가 CRM 옵션을 포함한 전체 지원 통합 목록은 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

## 설정

`sales/.claude/settings.local.json` 로컬 설정 파일을 만들어 개인화할 수 있습니다.

```json
{
  "name": "홍길동",
  "title": "Account Executive",
  "company": "회사명",
  "quota": {
    "annual": 1000000,
    "quarterly": 250000
  },
  "product": {
    "name": "제품명",
    "value_props": [
      "핵심 가치 제안 1",
      "핵심 가치 제안 2"
    ],
    "competitors": [
      "경쟁사 A",
      "경쟁사 B"
    ]
  }
}
```

설정되어 있지 않으면 플러그인이 대화형으로 정보를 요청합니다.
