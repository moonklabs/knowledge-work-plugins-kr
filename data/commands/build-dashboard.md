---
description: 차트, 필터, 테이블이 포함된 대화형 HTML 대시보드 구축
argument-hint: "<설명> [데이터 소스]"
---

# /build-dashboard - 대화형 대시보드 구축

> 익숙하지 않은 플레이스홀더가 있거나 연결된 도구를 확인해야 하는 경우 [CONNECTORS.md](../CONNECTORS.md)를 참조하십시오.

차트, 필터, 테이블, 전문적인 스타일링이 포함된 자체 완결형 대화형 HTML 대시보드를 구축합니다. 서버나 의존성 없이 브라우저에서 직접 열 수 있습니다.

## 사용법

```
/build-dashboard <대시보드 설명> [데이터 소스]
```

## 워크플로우

### 1. 대시보드 요구사항 이해

다음을 결정합니다:

- **목적**: 임원 개요, 운영 모니터링, 심층 분석, 팀 보고
- **대상**: 이 대시보드를 누가 사용할 것인가?
- **핵심 지표**: 가장 중요한 수치는 무엇인가?
- **차원**: 사용자가 어떤 기준으로 필터링하거나 분할할 수 있어야 하는가?
- **데이터 소스**: 실시간 쿼리, 붙여넣은 데이터, CSV 파일 또는 샘플 데이터

### 2. 데이터 수집

**데이터 웨어하우스가 연결된 경우:**
1. 필요한 데이터를 쿼리합니다
2. 결과를 HTML 파일 내에 JSON으로 포함합니다

**데이터가 붙여넣기 또는 업로드된 경우:**
1. 데이터를 파싱하고 정리합니다
2. 대시보드에 JSON으로 포함합니다

**데이터 없이 설명만으로 작업하는 경우:**
1. 설명된 스키마에 맞는 현실적인 샘플 데이터셋을 생성합니다
2. 대시보드에 샘플 데이터를 사용한다고 표시합니다
3. 실제 데이터로 교체하는 방법을 안내합니다

### 3. 대시보드 레이아웃 설계

표준 대시보드 레이아웃 패턴을 따릅니다:

```
┌──────────────────────────────────────────────────┐
│  Dashboard Title                    [Filters ▼]  │
├────────────┬────────────┬────────────┬───────────┤
│  KPI Card  │  KPI Card  │  KPI Card  │ KPI Card  │
├────────────┴────────────┼────────────┴───────────┤
│                         │                        │
│    Primary Chart        │   Secondary Chart      │
│    (largest area)       │                        │
│                         │                        │
├─────────────────────────┴────────────────────────┤
│                                                  │
│    Detail Table (sortable, scrollable)           │
│                                                  │
└──────────────────────────────────────────────────┘
```

**콘텐츠에 맞게 레이아웃을 조정합니다:**
- 상단에 2-4개의 KPI 카드로 헤드라인 수치 표시
- 중간 섹션에 1-3개의 차트로 추세와 분석 표시
- 하단에 드릴다운 데이터를 위한 선택적 상세 테이블
- 복잡도에 따라 헤더 또는 사이드바에 필터 배치

### 4. HTML 대시보드 구축

다음을 포함하는 단일 자체 완결형 HTML 파일을 생성합니다:

**구조 (HTML):**
- 시맨틱 HTML5 레이아웃
- CSS Grid 또는 Flexbox를 사용한 반응형 그리드
- 필터 컨트롤 (드롭다운, 날짜 선택기, 토글)
- 값과 레이블이 있는 KPI 카드
- 차트 컨테이너
- 정렬 가능한 헤더가 있는 데이터 테이블

**스타일링 (CSS):**
- 전문적인 색상 체계 (깔끔한 흰색, 회색, 데이터용 강조 색상)
- 은은한 그림자가 있는 카드 기반 레이아웃
- 일관된 타이포그래피 (빠른 로딩을 위한 시스템 폰트)
- 다양한 화면 크기에서 작동하는 반응형 디자인
- 인쇄 친화적 스타일

**상호작용 (JavaScript):**
- 대화형 차트를 위한 Chart.js (CDN을 통해 포함)
- 모든 차트와 테이블을 동시에 업데이트하는 필터 드롭다운
- 정렬 가능한 테이블 컬럼
- 차트의 호버 툴팁
- 숫자 형식 (쉼표, 통화, 백분율)

**데이터 (임베디드 JSON):**
- 모든 데이터를 HTML 내에 JavaScript 변수로 직접 포함
- 외부 데이터 페치 불필요
- 대시보드가 완전히 오프라인으로 작동

### 5. 차트 유형 구현

모든 차트에 Chart.js를 사용합니다. 일반적인 대시보드 차트 패턴:

- **꺾은선 차트**: 시계열 추세
- **막대 차트**: 카테고리 비교
- **도넛 차트**: 구성 비율 (6개 미만 카테고리일 때)
- **누적 막대 차트**: 시간에 따른 구성 변화
- **혼합형 (막대 + 꺾은선)**: 볼륨과 비율 오버레이

### 6. 상호작용 추가

**필터:**
```javascript
// 모든 필터는 중앙 필터 상태를 업데이트합니다
// 필터가 변경되면 차트와 테이블이 다시 렌더링됩니다
function applyFilters() {
    const filtered = data.filter(row => matchesFilters(row));
    updateKPIs(filtered);
    updateCharts(filtered);
    updateTable(filtered);
}
```

**테이블 정렬:**
- 컬럼 헤더를 클릭하여 오름차순/내림차순 정렬
- 현재 정렬 컬럼과 방향에 대한 시각적 표시

**툴팁:**
- 차트에서 호버 시 상세 값 표시
- KPI 카드에서 이전 기간과의 비교 표시

### 7. 저장 및 열기

1. 대시보드를 설명적인 이름의 HTML 파일로 저장합니다 (예: `sales_dashboard.html`)
2. 사용자의 기본 브라우저에서 엽니다
3. 올바르게 렌더링되는지 확인합니다
4. 데이터 업데이트 또는 커스터마이징 방법을 안내합니다

## 출력 템플릿

생성된 HTML 파일은 다음 구조를 따릅니다:

```html
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>[Dashboard Title]</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js@4.5.1" integrity="sha384-jb8JQMbMoBUzgWatfe6COACi2ljcDdZQ2OxczGA3bGNeWe+6DChMTBJemed7ZnvJ" crossorigin="anonymous"></script>
    <style>
        /* Professional dashboard CSS */
    </style>
</head>
<body>
    <div class="dashboard">
        <header><!-- Title and filters --></header>
        <section class="kpis"><!-- KPI cards --></section>
        <section class="charts"><!-- Chart containers --></section>
        <section class="details"><!-- Data table --></section>
    </div>
    <script>
        const DATA = [/* embedded JSON data */];
        // Dashboard initialization and interactivity
    </script>
</body>
</html>
```

## 예시

```
/build-dashboard 수익 추세, 주요 제품, 지역별 세부 내역을 포함한 월간 매출 대시보드를 만들어 주세요. 데이터는 orders 테이블에 있습니다.
```

```
/build-dashboard 여기 지원 티켓 데이터가 있습니다 [CSV 붙여넣기]. 우선순위별 양, 응답 시간 추세, 해결률을 보여주는 대시보드를 구축해 주세요.
```

```
/build-dashboard MRR, 이탈률, 신규 고객, NPS를 보여주는 SaaS 기업용 템플릿 임원 대시보드를 만들어 주세요. 샘플 데이터를 사용하세요.
```

## 팁

- 대시보드는 완전히 자체 완결형 HTML 파일입니다 -- 파일을 전송하여 누구에게나 공유할 수 있습니다
- 실시간 대시보드의 경우 BI 도구 연결을 고려하십시오. 이 대시보드는 특정 시점의 스냅샷입니다
- 다른 스타일링을 위해 "다크 모드" 또는 "프레젠테이션 모드"를 요청할 수 있습니다
- 브랜드에 맞는 특정 색상 체계를 요청할 수 있습니다
