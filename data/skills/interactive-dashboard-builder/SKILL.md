---
name: interactive-dashboard-builder
description: Chart.js, 드롭다운 필터, 전문적인 스타일링을 갖춘 자체 완결형 대화형 HTML 대시보드를 구축합니다. 대시보드를 생성하거나, 대화형 보고서를 구축하거나, 서버 없이 작동하는 차트와 필터가 포함된 공유 가능한 HTML 파일을 생성할 때 사용합니다.
---

# 대화형 대시보드 빌더 스킬

Chart.js, 필터, 상호작용, 전문적인 스타일링을 갖춘 자체 완결형 HTML/JS 대시보드 구축을 위한 패턴과 기법입니다.

## HTML/JS 대시보드 패턴

### 기본 템플릿

모든 대시보드는 다음 구조를 따릅니다:

```html
<!DOCTYPE html>
<html lang="ko">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>대시보드 제목</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js@4.5.1" integrity="sha384-jb8JQMbMoBUzgWatfe6COACi2ljcDdZQ2OxczGA3bGNeWe+6DChMTBJemed7ZnvJ" crossorigin="anonymous"></script>
    <script src="https://cdn.jsdelivr.net/npm/chartjs-adapter-date-fns@3.0.0" integrity="sha384-cVMg8E3QFwTvGCDuK+ET4PD341jF3W8nO1auiXfuZNQkzbUUiBGLsIQUE+b1mxws" crossorigin="anonymous"></script>
    <style>
        /* 대시보드 스타일 */
    </style>
</head>
<body>
    <div class="dashboard-container">
        <header class="dashboard-header">
            <h1>대시보드 제목</h1>
            <div class="filters">
                <!-- 필터 컨트롤 -->
            </div>
        </header>

        <section class="kpi-row">
            <!-- KPI 카드 -->
        </section>

        <section class="chart-row">
            <!-- 차트 컨테이너 -->
        </section>

        <section class="table-section">
            <!-- 데이터 테이블 -->
        </section>

        <footer class="dashboard-footer">
            <span>데이터 기준일: <span id="data-date"></span></span>
        </footer>
    </div>

    <script>
        // 포함된 데이터
        const DATA = [];

        // 대시보드 로직
        class Dashboard {
            constructor(data) {
                this.rawData = data;
                this.filteredData = data;
                this.charts = {};
                this.init();
            }

            init() {
                this.setupFilters();
                this.renderKPIs();
                this.renderCharts();
                this.renderTable();
            }

            applyFilters() {
                // 필터 로직
                this.filteredData = this.rawData.filter(row => {
                    // 각 활성 필터 적용
                    return true; // 플레이스홀더
                });
                this.renderKPIs();
                this.updateCharts();
                this.renderTable();
            }

            // ... 각 섹션별 메서드
        }

        const dashboard = new Dashboard(DATA);
    </script>
</body>
</html>
```

### KPI 카드 패턴

```html
<div class="kpi-card">
    <div class="kpi-label">총 매출</div>
    <div class="kpi-value" id="kpi-revenue">$0</div>
    <div class="kpi-change positive" id="kpi-revenue-change">+0%</div>
</div>
```

```javascript
function renderKPI(elementId, value, previousValue, format = 'number') {
    const el = document.getElementById(elementId);
    const changeEl = document.getElementById(elementId + '-change');

    // 값 형식 지정
    el.textContent = formatValue(value, format);

    // 증감 계산 및 표시
    if (previousValue && previousValue !== 0) {
        const pctChange = ((value - previousValue) / previousValue) * 100;
        const sign = pctChange >= 0 ? '+' : '';
        changeEl.textContent = `${sign}${pctChange.toFixed(1)}% vs 이전 기간`;
        changeEl.className = `kpi-change ${pctChange >= 0 ? 'positive' : 'negative'}`;
    }
}

function formatValue(value, format) {
    switch (format) {
        case 'currency':
            if (value >= 1e6) return `$${(value / 1e6).toFixed(1)}M`;
            if (value >= 1e3) return `$${(value / 1e3).toFixed(1)}K`;
            return `$${value.toFixed(0)}`;
        case 'percent':
            return `${value.toFixed(1)}%`;
        case 'number':
            if (value >= 1e6) return `${(value / 1e6).toFixed(1)}M`;
            if (value >= 1e3) return `${(value / 1e3).toFixed(1)}K`;
            return value.toLocaleString();
        default:
            return value.toString();
    }
}
```

### 차트 컨테이너 패턴

```html
<div class="chart-container">
    <h3 class="chart-title">월간 수익 추세</h3>
    <canvas id="revenue-chart"></canvas>
</div>
```

## Chart.js 통합

### 꺾은선 차트

```javascript
function createLineChart(canvasId, labels, datasets) {
    const ctx = document.getElementById(canvasId).getContext('2d');
    return new Chart(ctx, {
        type: 'line',
        data: {
            labels: labels,
            datasets: datasets.map((ds, i) => ({
                label: ds.label,
                data: ds.data,
                borderColor: COLORS[i % COLORS.length],
                backgroundColor: COLORS[i % COLORS.length] + '20',
                borderWidth: 2,
                fill: ds.fill || false,
                tension: 0.3,
                pointRadius: 3,
                pointHoverRadius: 6,
            }))
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            interaction: {
                mode: 'index',
                intersect: false,
            },
            plugins: {
                legend: {
                    position: 'top',
                    labels: { usePointStyle: true, padding: 20 }
                },
                tooltip: {
                    callbacks: {
                        label: function(context) {
                            return `${context.dataset.label}: ${formatValue(context.parsed.y, 'currency')}`;
                        }
                    }
                }
            },
            scales: {
                x: {
                    grid: { display: false }
                },
                y: {
                    beginAtZero: true,
                    ticks: {
                        callback: function(value) {
                            return formatValue(value, 'currency');
                        }
                    }
                }
            }
        }
    });
}
```

### 막대 차트

```javascript
function createBarChart(canvasId, labels, data, options = {}) {
    const ctx = document.getElementById(canvasId).getContext('2d');
    const isHorizontal = options.horizontal || labels.length > 8;

    return new Chart(ctx, {
        type: 'bar',
        data: {
            labels: labels,
            datasets: [{
                label: options.label || '값',
                data: data,
                backgroundColor: options.colors || COLORS.map(c => c + 'CC'),
                borderColor: options.colors || COLORS,
                borderWidth: 1,
                borderRadius: 4,
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            indexAxis: isHorizontal ? 'y' : 'x',
            plugins: {
                legend: { display: false },
                tooltip: {
                    callbacks: {
                        label: function(context) {
                            return formatValue(context.parsed[isHorizontal ? 'x' : 'y'], options.format || 'number');
                        }
                    }
                }
            },
            scales: {
                x: {
                    beginAtZero: true,
                    grid: { display: isHorizontal },
                    ticks: isHorizontal ? {
                        callback: function(value) {
                            return formatValue(value, options.format || 'number');
                        }
                    } : {}
                },
                y: {
                    beginAtZero: !isHorizontal,
                    grid: { display: !isHorizontal },
                    ticks: !isHorizontal ? {
                        callback: function(value) {
                            return formatValue(value, options.format || 'number');
                        }
                    } : {}
                }
            }
        }
    });
}
```

### 도넛 차트

```javascript
function createDoughnutChart(canvasId, labels, data) {
    const ctx = document.getElementById(canvasId).getContext('2d');
    return new Chart(ctx, {
        type: 'doughnut',
        data: {
            labels: labels,
            datasets: [{
                data: data,
                backgroundColor: COLORS.map(c => c + 'CC'),
                borderColor: '#ffffff',
                borderWidth: 2,
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            cutout: '60%',
            plugins: {
                legend: {
                    position: 'right',
                    labels: { usePointStyle: true, padding: 15 }
                },
                tooltip: {
                    callbacks: {
                        label: function(context) {
                            const total = context.dataset.data.reduce((a, b) => a + b, 0);
                            const pct = ((context.parsed / total) * 100).toFixed(1);
                            return `${context.label}: ${formatValue(context.parsed, 'number')} (${pct}%)`;
                        }
                    }
                }
            }
        }
    });
}
```

### 필터 변경 시 차트 업데이트

```javascript
function updateChart(chart, newLabels, newData) {
    chart.data.labels = newLabels;

    if (Array.isArray(newData[0])) {
        // 다중 데이터셋
        newData.forEach((data, i) => {
            chart.data.datasets[i].data = data;
        });
    } else {
        chart.data.datasets[0].data = newData;
    }

    chart.update('none'); // 'none'은 즉각적인 업데이트를 위해 애니메이션을 비활성화함
}
```

## 필터 및 상호작용 구현

### 드롭다운 필터

```html
<div class="filter-group">
    <label for="filter-region">지역</label>
    <select id="filter-region" onchange="dashboard.applyFilters()">
        <option value="all">모든 지역</option>
    </select>
</div>
```

```javascript
function populateFilter(selectId, data, field) {
    const select = document.getElementById(selectId);
    const values = [...new Set(data.map(d => d[field]))].sort();

    // "모두" 옵션 유지, 고유값 추가
    values.forEach(val => {
        const option = document.createElement('option');
        option.value = val;
        option.textContent = val;
        select.appendChild(option);
    });
}

function getFilterValue(selectId) {
    const val = document.getElementById(selectId).value;
    return val === 'all' ? null : val;
}
```

### 날짜 범위 필터

```html
<div class="filter-group">
    <label>날짜 범위</label>
    <input type="date" id="filter-date-start" onchange="dashboard.applyFilters()">
    <span>~</span>
    <input type="date" id="filter-date-end" onchange="dashboard.applyFilters()">
</div>
```

```javascript
function filterByDateRange(data, dateField, startDate, endDate) {
    return data.filter(row => {
        const rowDate = new Date(row[dateField]);
        if (startDate && rowDate < new Date(startDate)) return false;
        if (endDate && rowDate > new Date(endDate)) return false;
        return true;
    });
}
```

### 결합된 필터 로직

```javascript
applyFilters() {
    const region = getFilterValue('filter-region');
    const category = getFilterValue('filter-category');
    const startDate = document.getElementById('filter-date-start').value;
    const endDate = document.getElementById('filter-date-end').value;

    this.filteredData = this.rawData.filter(row => {
        if (region && row.region !== region) return false;
        if (category && row.category !== category) return false;
        if (startDate && row.date < startDate) return false;
        if (endDate && row.date > endDate) return false;
        return true;
    });

    this.renderKPIs();
    this.updateCharts();
    this.renderTable();
}
```

### 정렬 가능한 테이블

```javascript
function renderTable(containerId, data, columns) {
    const container = document.getElementById(containerId);
    let sortCol = null;
    let sortDir = 'desc';

    function render(sortedData) {
        let html = '<table class="data-table">';

        // 헤더
        html += '<thead><tr>';
        columns.forEach(col => {
            const arrow = sortCol === col.field
                ? (sortDir === 'asc' ? ' ▲' : ' ▼')
                : '';
            html += `<th onclick="sortTable('${col.field}')" style="cursor:pointer">${col.label}${arrow}</th>`;
        });
        html += '</tr></thead>';

        // 본문
        html += '<tbody>';
        sortedData.forEach(row => {
            html += '<tr>';
            columns.forEach(col => {
                const value = col.format ? formatValue(row[col.field], col.format) : row[col.field];
                html += `<td>${value}</td>`;
            });
            html += '</td>';
        });
        html += '</tbody></table>';

        container.innerHTML = html;
    }

    window.sortTable = function(field) {
        if (sortCol === field) {
            sortDir = sortDir === 'asc' ? 'desc' : 'asc';
        } else {
            sortCol = field;
            sortDir = 'desc';
        }
        const sorted = [...data].sort((a, b) => {
            const aVal = a[field], bVal = b[field];
            const cmp = aVal < bVal ? -1 : aVal > bVal ? 1 : 0;
            return sortDir === 'asc' ? cmp : -cmp;
        });
        render(sorted);
    };

    render(data);
}
```

## 대시보드용 CSS 스타일링

### 색상 시스템

```css
:root {
    /* 배경 레이어 */
    --bg-primary: #f8f9fa;
    --bg-card: #ffffff;
    --bg-header: #1a1a2e;

    /* 텍스트 */
    --text-primary: #212529;
    --text-secondary: #6c757d;
    --text-on-dark: #ffffff;

    /* 데이터용 강조 색상 */
    --color-1: #4C72B0;
    --color-2: #DD8452;
    --color-3: #55A868;
    --color-4: #C44E52;
    --color-5: #8172B3;
    --color-6: #937860;

    /* 상태 색상 */
    --positive: #28a745;
    --negative: #dc3545;
    --neutral: #6c757d;

    /* 간격 */
    --gap: 16px;
    --radius: 8px;
}
```

### 레이아웃

```css
* {
    margin: 0;
    padding: 0;
    box-sizing: border-box;
}

body {
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    background: var(--bg-primary);
    color: var(--text-primary);
    line-height: 1.5;
}

.dashboard-container {
    max-width: 1400px;
    margin: 0 auto;
    padding: var(--gap);
}

.dashboard-header {
    background: var(--bg-header);
    color: var(--text-on-dark);
    padding: 20px 24px;
    border-radius: var(--radius);
    margin-bottom: var(--gap);
    display: flex;
    justify-content: space-between;
    align-items: center;
    flex-wrap: wrap;
    gap: 12px;
}

.dashboard-header h1 {
    font-size: 20px;
    font-weight: 600;
}
```

### KPI 카드

```css
.kpi-row {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
    gap: var(--gap);
    margin-bottom: var(--gap);
}

.kpi-card {
    background: var(--bg-card);
    border-radius: var(--radius);
    padding: 20px 24px;
    box-shadow: 0 1px 3px rgba(0, 0, 0, 0.08);
}

.kpi-label {
    font-size: 13px;
    color: var(--text-secondary);
    text-transform: uppercase;
    letter-spacing: 0.5px;
    margin-bottom: 4px;
}

.kpi-value {
    font-size: 28px;
    font-weight: 700;
    color: var(--text-primary);
    margin-bottom: 4px;
}

.kpi-change {
    font-size: 13px;
    font-weight: 500;
}

.kpi-change.positive { color: var(--positive); }
.kpi-change.negative { color: var(--negative); }
```

### 차트 컨테이너

```css
.chart-row {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
    gap: var(--gap);
    margin-bottom: var(--gap);
}

.chart-container {
    background: var(--bg-card);
    border-radius: var(--radius);
    padding: 20px 24px;
    box-shadow: 0 1px 3px rgba(0, 0, 0, 0.08);
}

.chart-container h3 {
    font-size: 14px;
    font-weight: 600;
    color: var(--text-primary);
    margin-bottom: 16px;
}

.chart-container canvas {
    max-height: 300px;
}
```

### 필터

```css
.filters {
    display: flex;
    gap: 12px;
    align-items: center;
    flex-wrap: wrap;
}

.filter-group {
    display: flex;
    align-items: center;
    gap: 6px;
}

.filter-group label {
    font-size: 12px;
    color: rgba(255, 255, 255, 0.7);
}

.filter-group select,
.filter-group input[type="date"] {
    padding: 6px 10px;
    border: 1px solid rgba(255, 255, 255, 0.2);
    border-radius: 4px;
    background: rgba(255, 255, 255, 0.1);
    color: var(--text-on-dark);
    font-size: 13px;
}

.filter-group select option {
    background: var(--bg-header);
    color: var(--text-on-dark);
}
```

### 데이터 테이블

```css
.table-section {
    background: var(--bg-card);
    border-radius: var(--radius);
    padding: 20px 24px;
    box-shadow: 0 1px 3px rgba(0, 0, 0, 0.08);
    overflow-x: auto;
}

.data-table {
    width: 100%;
    border-collapse: collapse;
    font-size: 13px;
}

.data-table thead th {
    text-align: left;
    padding: 10px 12px;
    border-bottom: 2px solid #dee2e6;
    color: var(--text-secondary);
    font-weight: 600;
    font-size: 12px;
    text-transform: uppercase;
    letter-spacing: 0.5px;
    white-space: nowrap;
    user-select: none;
}

.data-table thead th:hover {
    color: var(--text-primary);
    background: #f8f9fa;
}

.data-table tbody td {
    padding: 10px 12px;
    border-bottom: 1px solid #f0f0f0;
}

.data-table tbody tr:hover {
    background: #f8f9fa;
}

.data-table tbody tr:last-child td {
    border-bottom: none;
}
```

### 반응형 디자인

```css
@media (max-width: 768px) {
    .dashboard-header {
        flex-direction: column;
        align-items: flex-start;
    }

    .kpi-row {
        grid-template-columns: repeat(2, 1fr);
    }

    .chart-row {
        grid-template-columns: 1fr;
    }

    .filters {
        flex-direction: column;
        align-items: flex-start;
    }
}

@media print {
    body { background: white; }
    .dashboard-container { max-width: none; }
    .filters { display: none; }
    .chart-container { break-inside: avoid; }
    .kpi-card { border: 1px solid #dee2e6; box-shadow: none; }
}
```

## 대용량 데이터셋의 성능 고려사항

### 데이터 크기 가이드라인

| 데이터 크기 | 접근 방법 |
|---|---|
| 1,000행 미만 | HTML에 직접 포함. 완전한 상호작용. |
| 1,000 - 10,000행 | HTML에 포함. 차트용으로 사전 집계가 필요할 수 있음. |
| 10,000 - 100,000행 | 서버 측에서 사전 집계. 집계된 데이터만 포함. |
| 100,000행 초과 | 클라이언트 측 대시보드에 부적합. BI 도구를 사용하거나 페이지네이션. |

### 사전 집계 패턴

원시 데이터를 포함하고 브라우저에서 집계하는 대신:

```javascript
// X: 50,000개의 원시 행을 포함하지 마십시오
const RAW_DATA = [/* 50,000 rows */];

// O: 포함하기 전에 결과 데이터를 사전 집계하십시오
const CHART_DATA = {
    monthly_revenue: [
        { month: '2024-01', revenue: 150000, orders: 1200 },
        { month: '2024-02', revenue: 165000, orders: 1350 },
        // ... 50,000행 대신 12행
    ],
    top_products: [
        { product: 'Widget A', revenue: 45000 },
        // ... 10행
    ],
    kpis: {
        total_revenue: 1980000,
        total_orders: 15600,
        avg_order_value: 127,
    }
};
```

### 차트 성능

- 꺾은선 차트는 시리즈당 500개 미만의 데이터 포인트로 제한합니다 (필요시 다운샘플링)
- 막대 차트는 50개 미만의 카테고리로 제한합니다
- 산점도는 1,000개 포인트로 제한합니다 (대규모 데이터셋에는 샘플링 사용)
- 차트가 많은 대시보드에서는 애니메이션을 비활성화합니다: Chart.js 옵션에서 `animation: false`
- 필터 트리거 업데이트에는 `Chart.update()` 대신 `Chart.update('none')`을 사용합니다

### DOM 성능

- 데이터 테이블은 100-200개 가시 행으로 제한합니다. 더 많은 경우 페이지네이션을 추가합니다.
- 조율된 차트 업데이트에는 `requestAnimationFrame`을 사용합니다
- 필터 변경 시 전체 DOM을 다시 구축하지 않습니다 -- 변경된 요소만 업데이트합니다

```javascript
// 효율적인 테이블 페이지네이션
function renderTablePage(data, page, pageSize = 50) {
    const start = page * pageSize;
    const end = Math.min(start + pageSize, data.length);
    const pageData = data.slice(start, end);
    // pageData만 렌더링
    // 페이지네이션 컨트롤 표시: "2,340개 중 1-50개 표시 중"
}
```
