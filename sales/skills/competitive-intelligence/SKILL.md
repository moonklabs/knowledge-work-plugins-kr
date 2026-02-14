---
name: competitive-intelligence
description: 경쟁사를 조사하고 인터랙티브 배틀카드를 구축합니다. 클릭 가능한 경쟁사 카드와 비교 매트릭스가 포함된 HTML 아티팩트를 출력합니다. "competitive intel", "research competitors", "how do we compare to [경쟁사]", "battlecard for [경쟁사]", "what's new with [경쟁사]" 등으로 실행합니다.
---

# 경쟁 정보

경쟁사를 광범위하게 조사하고 거래에서 활용할 수 있는 **인터랙티브 HTML 배틀카드**를 생성합니다. 출력물은 클릭 가능한 경쟁사 탭과 전체 비교 매트릭스가 포함된 독립형 아티팩트입니다.

## 작동 방식

```
┌─────────────────────────────────────────────────────────────────┐
│                  COMPETITIVE INTELLIGENCE                        │
├─────────────────────────────────────────────────────────────────┤
│  기본 기능 (웹 검색으로 단독 작동)                                  │
│  ✓ 경쟁사 제품 심층 분석: 기능, 가격, 포지셔닝                      │
│  ✓ 최근 출시: 최근 90일간 출시한 제품                              │
│  ✓ 자사 출시: 대응으로 출시한 제품                                 │
│  ✓ 차별화 매트릭스: 우리가 이기는 곳 vs. 지는 곳                    │
│  ✓ 세일즈 토크 트랙: 각 경쟁사 대비 포지셔닝 방법                   │
│  ✓ 랜드마인 질문: 약점을 자연스럽게 노출시키는 질문                  │
├─────────────────────────────────────────────────────────────────┤
│  출력: 인터랙티브 HTML 배틀카드                                    │
│  ✓ 비교 매트릭스 개요                                             │
│  ✓ 각 경쟁사별 클릭 가능한 탭                                     │
│  ✓ 다크 테마, 전문적 스타일링                                     │
│  ✓ 독립형 HTML 파일 — 어디서든 공유 또는 호스팅 가능                │
├─────────────────────────────────────────────────────────────────┤
│  강화 모드 (도구 연결 시)                                         │
│  + CRM: 성패 데이터, 종료된 거래의 경쟁사 언급                      │
│  + 문서: 기존 배틀카드, 경쟁 전략 플레이북                          │
│  + 채팅: 내부 인텔리전스, 동료의 현장 보고                          │
│  + 녹취록: 고객 통화 중 경쟁사 언급                                │
└─────────────────────────────────────────────────────────────────┘
```

---

## 시작하기

이 스킬 실행 시 맥락을 요청합니다:

**필수:**
- 어떤 기업에서 일하시나요? (또는 이메일에서 자동 감지)
- 주요 경쟁사는 누구인가요? (1-5개)

**선택사항:**
- 어떤 경쟁사를 먼저 집중적으로 분석할까요?
- 경쟁 중인 특정 거래가 있나요?
- 경쟁사에 대해 고객으로부터 들은 페인 포인트가 있나요?

이전 세션에서 판매자 맥락이 있으면 확인 후 질문을 건너뜁니다.

---

## 커넥터 (선택사항)

| 커넥터 | 추가 기능 |
|--------|----------|
| **CRM** | 각 경쟁사별 성패 이력, 거래 수준 경쟁사 추적 |
| **문서** | 기존 배틀카드, 제품 비교 문서, 경쟁 전략 플레이북 |
| **채팅** | 내부 채팅 인텔리전스(예: Slack) — 팀이 현장에서 듣는 정보 |
| **녹취록** | 고객 통화 중 경쟁사 언급, 제기된 이의 |

> **커넥터가 없나요?** 웹 리서치만으로도 훌륭합니다. 제품 페이지, 가격, 블로그, 릴리스 노트, 리뷰, 채용 공고 등 공개 소스에서 모든 정보를 수집합니다.

---

## 출력: 인터랙티브 HTML 배틀카드

이 스킬은 다음을 포함하는 **독립형 HTML 파일**을 생성합니다:

### 1. 비교 매트릭스 (메인 뷰)
자사 vs. 모든 경쟁사를 한눈에 비교하는 개요:
- 기능 비교 그리드
- 가격 비교
- 시장 포지셔닝
- 성사율 지표 (CRM 연결 시)

### 2. 경쟁사 탭 (클릭하여 확장)
각 경쟁사별 클릭 가능한 카드로 다음을 표시:
- 기업 프로필 (규모, 투자, 타겟 시장)
- 판매 제품 및 포지셔닝 방식
- 최근 출시 (최근 90일)
- 경쟁사가 이기는 곳 vs. 우리가 이기는 곳
- 가격 인텔리전스
- 시나리오별 토크 트랙
- 이의 처리
- 랜드마인 질문

### 3. 자사 카드
- 자사 출시 (최근 90일)
- 핵심 차별화 요소
- 증거 및 고객 인용

---

## HTML 구조

```html
<!DOCTYPE html>
<html>
<head>
    <title>Battlecard: [Your Company] vs Competitors</title>
    <style>
        /* Dark theme, professional styling */
        /* Tabbed navigation */
        /* Expandable cards */
        /* Responsive design */
    </style>
</head>
<body>
    <!-- Header with your company + date -->
    <header>
        <h1>[Your Company] Competitive Battlecard</h1>
        <p>Generated: [Date] | Competitors: [List]</p>
    </header>

    <!-- Tab Navigation -->
    <nav class="tabs">
        <button class="tab active" data-tab="matrix">Comparison Matrix</button>
        <button class="tab" data-tab="competitor-1">[Competitor 1]</button>
        <button class="tab" data-tab="competitor-2">[Competitor 2]</button>
        <button class="tab" data-tab="competitor-3">[Competitor 3]</button>
    </nav>

    <!-- Comparison Matrix Tab -->
    <section id="matrix" class="tab-content active">
        <h2>Head-to-Head Comparison</h2>
        <table class="comparison-matrix">
            <!-- Feature rows with you vs each competitor -->
        </table>

        <h2>Quick Win/Loss Guide</h2>
        <div class="win-loss-grid">
            <!-- Per-competitor: when you win, when you lose -->
        </div>
    </section>

    <!-- Individual Competitor Tabs -->
    <section id="competitor-1" class="tab-content">
        <div class="battlecard">
            <div class="profile"><!-- Company info --></div>
            <div class="differentiation"><!-- Where they win / you win --></div>
            <div class="talk-tracks"><!-- Scenario-based positioning --></div>
            <div class="objections"><!-- Common objections + responses --></div>
            <div class="landmines"><!-- Questions to ask --></div>
        </div>
    </section>

    <script>
        // Tab switching logic
        // Expand/collapse sections
    </script>
</body>
</html>
```

---

## 비주얼 디자인

### 컬러 시스템
```css
:root {
    /* Dark theme base */
    --bg-primary: #0a0d14;
    --bg-elevated: #0f131c;
    --bg-surface: #161b28;
    --bg-hover: #1e2536;

    /* Text */
    --text-primary: #ffffff;
    --text-secondary: rgba(255, 255, 255, 0.7);
    --text-muted: rgba(255, 255, 255, 0.5);

    /* Accent (your brand or neutral) */
    --accent: #3b82f6;
    --accent-hover: #2563eb;

    /* Status indicators */
    --you-win: #10b981;
    --they-win: #ef4444;
    --tie: #f59e0b;
}
```

### 카드 디자인
- 둥근 모서리 (12px)
- 미세한 테두리 (1px, 낮은 불투명도)
- 약간의 부양 효과가 있는 호버 상태
- 부드러운 전환 (200ms)

### 비교 매트릭스
- 고정 헤더 행
- 색상 코드 우승 지표 (녹색 = 자사, 빨간색 = 경쟁사, 노란색 = 무승부)
- 세부 정보 확장 가능한 행

---

## 실행 흐름

### 1단계: 판매자 맥락 수집

```
첫 사용 시:
1. 질문: "어떤 기업에서 일하시나요?"
2. 질문: "무엇을 판매하시나요? (제품/서비스 한 줄)"
3. 질문: "주요 경쟁사는 누구인가요? (최대 5개)"
4. 향후 세션을 위해 맥락 저장

재방문 사용자:
1. 확인: "아직 [기업]에서 [제품]을 판매하시나요?"
2. 질문: "동일한 경쟁사인가요, 추가할 새 경쟁사가 있나요?"
```

### 2단계: 자사 리서치 (항상 실행)

```
웹 검색:
1. "[자사] product" — 현재 제공 제품
2. "[자사] pricing" — 가격 모델
3. "[자사] news" — 최근 발표 (90일)
4. "[자사] product updates OR changelog OR releases" — 출시 현황
5. "[자사] vs [경쟁사]" — 기존 비교
```

### 3단계: 각 경쟁사 리서치 (항상 실행)

```
각 경쟁사에 대해 실행:
1. "[경쟁사] product features" — 제공 기능
2. "[경쟁사] pricing" — 과금 방식
3. "[경쟁사] news" — 최근 발표
4. "[경쟁사] product updates OR changelog OR releases" — 출시 현황
5. "[경쟁사] reviews G2 OR Capterra OR TrustRadius" — 고객 평가
6. "[경쟁사] vs [대안]" — 포지셔닝 방식
7. "[경쟁사] customers" — 고객층
8. "[경쟁사] careers" — 채용 신호 (성장 영역)
```

### 4단계: 연결된 소스 가져오기 (사용 가능 시)

```
CRM 연결 시:
1. 경쟁사 필드 = [경쟁사]인 성사 거래 조회
2. 경쟁사 필드 = [경쟁사]인 실패 거래 조회
3. 성패 패턴 추출

문서 연결 시:
1. "battlecard [경쟁사]" 검색
2. "competitive [경쟁사]" 검색
3. 기존 포지셔닝 문서 가져오기

채팅 연결 시:
1. "[경쟁사]" 언급 검색 (최근 90일)
2. 현장 인텔리전스 및 동료 인사이트 추출

녹취록 연결 시:
1. 통화에서 "[경쟁사]" 언급 검색
2. 이의 및 고객 인용 추출
```

### 5단계: HTML 아티팩트 구축

```
1. 각 경쟁사별 데이터 구조화
2. 비교 매트릭스 구축
3. 개별 배틀카드 생성
4. 시나리오별 토크 트랙 생성
5. 랜드마인 질문 정리
6. 독립형 HTML로 렌더링
7. [자사]-battlecard-[날짜].html로 저장
```

---

## 경쟁사별 데이터 구조

```yaml
competitor:
  name: "[이름]"
  website: "[URL]"
  profile:
    founded: "[연도]"
    funding: "[단계 + 금액]"
    employees: "[인원]"
    target_market: "[판매 대상]"
    pricing_model: "[인당 / 사용량 기반 / 등]"
    market_position: "[리더 / 챌린저 / 니치]"

  what_they_sell: "[제품 요약]"
  their_positioning: "[자체 설명 방식]"

  recent_releases:
    - date: "[날짜]"
      release: "[기능/제품]"
      impact: "[중요한 이유]"

  where_they_win:
    - area: "[영역]"
      advantage: "[그들의 강점]"
      how_to_handle: "[대응 방법]"

  where_you_win:
    - area: "[영역]"
      advantage: "[우리의 강점]"
      proof_point: "[증거]"

  pricing:
    model: "[과금 방식]"
    entry_price: "[시작 가격]"
    enterprise: "[엔터프라이즈 가격]"
    hidden_costs: "[구현 비용 등]"
    talk_track: "[가격 논의 방법]"

  talk_tracks:
    early_mention: "[초기 언급 시 전략]"
    displacement: "[고객이 사용 중일 때 전략]"
    late_addition: "[평가 후반에 추가될 때 전략]"

  objections:
    - objection: "[고객 발언]"
      response: "[대응 방법]"

  landmines:
    - "[약점을 노출시키는 질문]"

  win_loss: # CRM 연결 시
    win_rate: "[X]%"
    common_win_factors: "[성공 예측 요인]"
    common_loss_factors: "[실패 예측 요인]"
```

---

## 전달

```markdown
## 배틀카드 생성 완료

[배틀카드 보기](file:///path/to/[YourCompany]-battlecard-[date].html)

---

**요약**
- **자사**: [이름]
- **분석 경쟁사**: [목록]
- **데이터 소스**: 웹 리서치 [+ CRM] [+ 문서] [+ 녹취록]

---

**활용 방법**
- **통화 전**: 관련 경쟁사 탭을 열고 토크 트랙 검토
- **통화 중**: 랜드마인 질문 참조
- **성패 후**: 새로운 인텔리전스로 업데이트

---

**공유 옵션**
- **로컬 파일**: 모든 브라우저에서 열기
- **호스팅**: Netlify, Vercel 또는 내부 위키에 업로드
- **직접 공유**: HTML 파일을 팀원에게 전송

---

**최신 상태 유지**
이 스킬을 다시 실행하면 최신 인텔리전스로 갱신됩니다. 권장: 월 1회 또는 주요 거래 전.
```

---

## 갱신 주기

경쟁 정보는 시간이 지나면 오래됩니다. 권장 갱신 주기:

| 트리거 | 조치 |
|--------|------|
| **월간** | 빠른 갱신 — 신규 출시, 뉴스, 가격 변경 |
| **주요 거래 전** | 해당 거래의 특정 경쟁사에 대한 심층 갱신 |
| **성패 후** | 새로운 데이터로 패턴 업데이트 |
| **경쟁사 발표 시** | 해당 경쟁사에 대한 즉시 업데이트 |

---

## 더 나은 인텔리전스를 위한 팁

1. **약점에 솔직하게** — 경쟁사의 강점을 인정하는 것이 신뢰를 만듭니다
2. **기능이 아닌 결과에 집중** — "X 기능이 있다"보다 "고객이 Y 결과를 달성한다"가 중요합니다
3. **현장에서 업데이트** — 최고의 인텔리전스는 웹사이트가 아닌 실제 고객 대화에서 나옵니다
4. **랜드마인 심기, 비방 금지** — 약점을 노출시키는 질문을 하되, 절대 비방하지 마세요
5. **출시 추적 철저히** — 경쟁사가 출시하는 것이 그들의 전략과 우리의 기회를 말해줍니다

---

## 관련 스킬

- **account-research** — 접촉 전 특정 잠재 고객 리서치
- **call-prep** — 경쟁사가 관련된 통화 준비
- **create-an-asset** — 특정 거래를 위한 맞춤형 비교 페이지 구축
