---
description: 차트, 필터, 테이블이 포함된 인터랙티브 HTML 대시보드 생성
argument-hint: "<description> [data source]"
---

# /build-dashboard - 인터랙티브 대시보드 생성

> 낯선 플레이스홀더가 보이거나 어떤 도구가 연결되어 있는지 확인하려면 [CONNECTORS.md](../CONNECTORS.md)를 보세요.

차트, 필터, 테이블, 전문적인 스타일을 포함한 자체 포함형 HTML 대시보드를 생성합니다. 브라우저에서 바로 열 수 있으며 서버나 의존성이 필요하지 않습니다.

## 사용법

```
/build-dashboard <대시보드 설명> [data source]
```

## 워크플로

### 1. 대시보드 요구사항 이해

다음을 확인합니다.

- **목적**: 임원 요약, 운영 모니터링, 심층 분석, 팀 리포팅
- **대상**: 누가 이 대시보드를 사용할 것인가
- **핵심 지표**: 가장 중요한 숫자
- **차원**: 사용자가 필터/분해할 차원
- **데이터 소스**: 라이브 쿼리, 붙여넣은 데이터, CSV 파일, 샘플 데이터

### 2. 데이터 수집

**데이터 웨어하우스가 연결된 경우:**
1. 필요한 데이터를 쿼리
2. 결과를 HTML 내부 JSON으로 임베드

**데이터가 붙여넣기/업로드된 경우:**
1. 데이터를 파싱하고 정제
2. JSON으로 대시보드에 임베드

**데이터 없이 설명만 있는 경우:**
1. 설명된 스키마에 맞는 현실적인 샘플 데이터 생성
2. 대시보드에 샘플 데이터 사용 사실을 명시
3. 실제 데이터로 교체하는 방법 안내

### 3. 대시보드 레이아웃 설계

표준 레이아웃 패턴을 따릅니다.

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

**콘텐츠에 맞게 레이아웃을 조정합니다.**
- 상단에 2~4개의 KPI 카드
- 중간 섹션에 1~3개의 차트(추세/분해)
- 하단에 드릴다운용 상세 테이블(선택)
- 필터는 복잡도에 따라 헤더 또는 사이드바

### 4. HTML 대시보드 구성

하나의 자체 포함형 HTML 파일을 생성합니다.

**구조(HTML):**
- 의미 있는 HTML5 레이아웃
- CSS Grid/Flexbox 기반 반응형 그리드
- 필터 컨트롤(드롭다운, 날짜 선택기, 토글)
- KPI 카드(값 + 레이블)
- 차트 컨테이너
- 정렬 가능한 데이터 테이블

**스타일(CSS):**
- 전문적인 컬러 스킴(깨끗한 흰색/회색 + 데이터 강조색)
- 카드 기반 레이아웃과 부드러운 그림자
- 일관된 타이포그래피(빠른 로딩을 위한 시스템 폰트)
- 다양한 화면 크기 대응
- 인쇄 친화 스타일

**인터랙션(JavaScript):**
- Chart.js(인터랙티브 차트, CDN 포함)
- 필터 변경 시 모든 차트/테이블 동시 갱신
- 테이블 컬럼 정렬
- 차트 툴팁
- 숫자 포맷(콤마, 통화, 퍼센트)

**데이터(임베디드 JSON):**
- 모든 데이터를 HTML에 직접 임베드
- 외부 데이터 fetch 없음
- 오프라인에서도 동작

### 5. 차트 유형 구현

모든 차트는 Chart.js를 사용합니다. 일반적인 패턴은 다음과 같습니다.

- **라인 차트**: 시간 추세
- **막대 차트**: 카테고리 비교
- **도넛 차트**: 구성(카테고리 <6일 때)
- **누적 막대**: 시간에 따른 구성
- **혼합형(막대 + 라인)**: 볼륨과 비율 동시 표시

### 6. 인터랙션 추가

**필터:**
```javascript
// 모든 필터는 중앙 필터 상태를 업데이트합니다
// 필터 변경 시 차트와 테이블을 다시 렌더링합니다
function applyFilters() {
    const filtered = data.filter(row => matchesFilters(row));
    updateKPIs(filtered);
    updateCharts(filtered);
    updateTable(filtered);
}
```

**테이블 정렬:**
- 컬럼 헤더 클릭 시 오름/내림차순 정렬
- 현재 정렬 컬럼과 방향 표시

**툴팁:**
- 차트 hover 시 상세 값 표시
- KPI 카드에 이전 기간 대비 표시

### 7. 저장 및 열기

1. 설명적인 이름으로 HTML 저장(예: `sales_dashboard.html`)
2. 기본 브라우저에서 열기
3. 렌더링 확인
4. 데이터 업데이트/커스터마이징 방법 안내

## 출력 템플릿

생성된 HTML 파일 구조는 다음과 같습니다.

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
/build-dashboard 월간 매출 추세, 상위 제품, 지역별 분해를 포함한 세일즈 대시보드. 데이터는 orders 테이블에 있습니다.
```

```
/build-dashboard 지원 티켓 데이터입니다 [CSV 붙여넣기]. 우선순위별 볼륨, 응답 시간 추세, 해결률을 보여주는 대시보드를 만들어 주세요.
```

```
/build-dashboard SaaS 회사를 위한 임원용 대시보드 템플릿을 만들어 주세요. MRR, 이탈, 신규 고객, NPS를 포함하고 샘플 데이터를 사용하세요.
```

## 팁

- 대시보드는 완전히 자체 포함된 HTML 파일입니다. 파일을 공유하면 누구나 볼 수 있습니다
- 실시간 대시보드가 필요하다면 BI 도구 연결을 고려하세요. 이 대시보드는 특정 시점 스냅샷입니다
- "dark mode" 또는 "presentation mode"를 요청하면 스타일을 조정합니다
- 브랜드에 맞는 특정 색상 스킴을 요청할 수 있습니다
