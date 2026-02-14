---
name: data-visualization
description: Python(matplotlib, seaborn, plotly)으로 효과적인 데이터 시각화를 생성합니다. 차트를 구축하거나, 데이터셋에 적합한 차트 유형을 선택하거나, 출판 수준의 그림을 생성하거나, 접근성 및 색상 이론과 같은 디자인 원칙을 적용할 때 사용합니다.
---

# 데이터 시각화 스킬

차트 선택 가이드, Python 시각화 코드 패턴, 디자인 원칙, 효과적인 데이터 시각화를 위한 접근성 고려사항입니다.

## 차트 선택 가이드

### 데이터 관계별 선택

| 표현하려는 내용 | 최적 차트 | 대안 |
|---|---|---|
| **시간에 따른 추세** | 꺾은선 차트 | 영역 차트 (누적 또는 구성 표시 시) |
| **카테고리 간 비교** | 세로 막대 차트 | 가로 막대 (카테고리가 많은 경우), 롤리팝 차트 |
| **순위** | 가로 막대 차트 | 점 그래프, 기울기 차트 (두 기간 비교) |
| **전체 대비 부분 구성** | 누적 막대 차트 | 트리맵 (계층적), 와플 차트 |
| **시간에 따른 구성** | 누적 영역 차트 | 100% 누적 막대 (비율 중심) |
| **분포** | 히스토그램 | 상자 그림 (그룹 비교), 바이올린 플롯, 스트립 플롯 |
| **상관관계 (2변수)** | 산점도 | 버블 차트 (3번째 변수를 크기로 추가) |
| **상관관계 (다변수)** | 히트맵 (상관 행렬) | 페어 플롯 |
| **지리적 패턴** | 코로플레스 맵 | 버블 맵, 헥스 맵 |
| **흐름/프로세스** | 생키 다이어그램 | 퍼널 차트 (순차적 단계) |
| **관계 네트워크** | 네트워크 그래프 | 코드 다이어그램 |
| **목표 대비 성과** | 불릿 차트 | 게이지 (단일 KPI만) |
| **여러 KPI 동시 표시** | 소규모 다중 차트 | 별도 차트가 있는 대시보드 |

### 특정 차트를 사용하지 말아야 할 때

- **파이 차트**: 6개 미만 카테고리이고 정확한 비율보다 대략적 비교가 중요한 경우가 아니면 지양. 인간은 각도 비교에 약합니다. 막대 차트를 대신 사용하십시오.
- **3D 차트**: 절대 사용하지 않습니다. 인식을 왜곡하고 정보를 추가하지 않습니다.
- **이중 축 차트**: 주의하여 사용합니다. 상관관계를 암시하여 오해를 줄 수 있습니다. 사용 시 양쪽 축에 명확한 레이블을 부여합니다.
- **누적 막대 (카테고리 다수)**: 중간 세그먼트 비교가 어렵습니다. 소규모 다중 또는 그룹 막대를 대신 사용합니다.
- **도넛 차트**: 파이 차트보다 약간 낫지만 근본적으로 같은 문제. 최대한 단일 KPI 표시용으로만 사용합니다.

## Python 시각화 코드 패턴

### 설정 및 스타일

```python
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
import pandas as pd
import numpy as np

# 전문 스타일 설정
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({
    'figure.figsize': (10, 6),
    'figure.dpi': 150,
    'font.size': 11,
    'axes.titlesize': 14,
    'axes.titleweight': 'bold',
    'axes.labelsize': 11,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 16,
})

# 색맹 친화적 팔레트
PALETTE_CATEGORICAL = ['#4C72B0', '#DD8452', '#55A868', '#C44E52', '#8172B3', '#937860']
PALETTE_SEQUENTIAL = 'YlOrRd'
PALETTE_DIVERGING = 'RdBu_r'
```

### 꺾은선 차트 (시계열)

```python
fig, ax = plt.subplots(figsize=(10, 6))

for label, group in df.groupby('category'):
    ax.plot(group['date'], group['value'], label=label, linewidth=2)

ax.set_title('카테고리별 지표 추세', fontweight='bold')
ax.set_xlabel('날짜')
ax.set_ylabel('값')
ax.legend(loc='upper left', frameon=True)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# x축 날짜 형식 지정
fig.autofmt_xdate()

plt.tight_layout()
plt.savefig('trend_chart.png', dpi=150, bbox_inches='tight')
```

### 막대 차트 (비교)

```python
fig, ax = plt.subplots(figsize=(10, 6))

# 읽기 쉽게 값 기준으로 정렬
df_sorted = df.sort_values('metric', ascending=True)

bars = ax.barh(df_sorted['category'], df_sorted['metric'], color=PALETTE_CATEGORICAL[0])

# 값 레이블 추가
for bar in bars:
    width = bar.get_width()
    ax.text(width + 0.5, bar.get_y() + bar.get_height()/2,
            f'{width:,.0f}', ha='left', va='center', fontsize=10)

ax.set_title('카테고리별 지표 (순위)', fontweight='bold')
ax.set_xlabel('지표 값')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('bar_chart.png', dpi=150, bbox_inches='tight')
```

### 히스토그램 (분포)

```python
fig, ax = plt.subplots(figsize=(10, 6))

ax.hist(df['value'], bins=30, color=PALETTE_CATEGORICAL[0], edgecolor='white', alpha=0.8)

# 평균 및 중앙값 선 추가
mean_val = df['value'].mean()
median_val = df['value'].median()
ax.axvline(mean_val, color='red', linestyle='--', linewidth=1.5, label=f'평균: {mean_val:,.1f}')
ax.axvline(median_val, color='green', linestyle='--', linewidth=1.5, label=f'중앙값: {median_val:,.1f}')

ax.set_title('값의 분포', fontweight='bold')
ax.set_xlabel('값')
ax.set_ylabel('빈도')
ax.legend()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('histogram.png', dpi=150, bbox_inches='tight')
```

### 히트맵

```python
fig, ax = plt.subplots(figsize=(10, 8))

# 히트맵 형식을 위한 데이터 피벗
pivot = df.pivot_table(index='row_dim', columns='col_dim', values='metric', aggfunc='sum')

sns.heatmap(pivot, annot=True, fmt=',.0f', cmap='YlOrRd',
            linewidths=0.5, ax=ax, cbar_kws={'label': '지표 값'})

ax.set_title('행 차원 및 열 차원 기준 지표', fontweight='bold')
ax.set_xlabel('열 차원')
ax.set_ylabel('행 차원')

plt.tight_layout()
plt.savefig('heatmap.png', dpi=150, bbox_inches='tight')
```

### 소규모 다중 차트

```python
categories = df['category'].unique()
n_cats = len(categories)
n_cols = min(3, n_cats)
n_rows = (n_cats + n_cols - 1) // n_cols

fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows), sharex=True, sharey=True)
axes = axes.flatten() if n_cats > 1 else [axes]

for i, cat in enumerate(categories):
    ax = axes[i]
    subset = df[df['category'] == cat]
    ax.plot(subset['date'], subset['value'], color=PALETTE_CATEGORICAL[i % len(PALETTE_CATEGORICAL)])
    ax.set_title(cat, fontsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# 빈 서브플롯 숨기기
for j in range(i+1, len(axes)):
    axes[j].set_visible(False)

fig.suptitle('카테고리별 추세', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('small_multiples.png', dpi=150, bbox_inches='tight')
```

### 숫자 형식 헬퍼

```python
def format_number(val, format_type='number'):
    """차트 레이블용 숫자 형식 지정"""
    if format_type == 'currency':
        if abs(val) >= 1e9:
            return f'${val/1e9:.1f}B'
        elif abs(val) >= 1e6:
            return f'${val/1e6:.1f}M'
        elif abs(val) >= 1e3:
            return f'${val/1e3:.1f}K'
        else:
            return f'${val:,.0f}'
    elif format_type == 'percent':
        return f'{val:.1f}%'
    elif format_type == 'number':
        if abs(val) >= 1e9:
            return f'{val/1e9:.1f}B'
        elif abs(val) >= 1e6:
            return f'{val/1e6:.1f}M'
        elif abs(val) >= 1e3:
            return f'{val/1e3:.1f}K'
        else:
            return f'{val:,.0f}'
    return str(val)

# 축 포맷터와 함께 사용
ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, p: format_number(x, 'currency')))
```

### Plotly를 활용한 대화형 차트

```python
import plotly.express as px
import plotly.graph_objects as go

# 간단한 대화형 꺾은선 차트
fig = px.line(df, x='date', y='value', color='category',
              title='대화형 지표 추세',
              labels={'value': '지표 값', 'date': '날짜'})
fig.update_layout(hovermode='x unified')
fig.write_html('interactive_chart.html')
fig.show()

# 상세 데이터가 포함된 대화형 산점도
fig = px.scatter(df, x='metric_a', y='metric_b', color='category',
                 size='size_metric', hover_data=['name', 'detail_field'],
                 title='상관관계 분석')
fig.show()
```

## 디자인 원칙

### 색상

- **색상을 목적 있게 사용**: 색상은 장식이 아닌 데이터를 인코딩해야 합니다
- **스토리를 강조**: 핵심 인사이트에 밝은 강조 색상을 사용하고, 나머지는 회색 처리합니다
- **순차적 데이터**: 정렬된 값에는 단일 색조 그라데이션(밝은색에서 어두운색)을 사용합니다
- **발산적 데이터**: 의미 있는 중심이 있는 데이터에는 중립 중간점이 있는 이중 색조 그라데이션을 사용합니다
- **범주형 데이터**: 구분되는 색조를 사용하되, 6-8개를 초과하면 혼란스러워집니다
- **적색/녹색만 사용 지양**: 남성의 8%가 적녹 색맹입니다. 주 색상 쌍으로 파란색/주황색을 사용합니다

### 타이포그래피

- **제목은 인사이트를 명시**: "매출이 전년 대비 23% 성장"이 "월별 매출"보다 좋습니다
- **부제는 맥락을 추가**: 날짜 범위, 적용된 필터, 데이터 소스
- **축 레이블은 읽기 쉽게**: 가능하면 90도 회전 지양. 단축하거나 줄바꿈합니다
- **데이터 레이블은 정밀도를 추가**: 핵심 포인트에만, 모든 막대에 달지 않습니다
- **주석으로 강조**: 텍스트 주석으로 특정 포인트를 호출합니다

### 레이아웃

- **차트 잡음을 줄입니다**: 정보를 전달하지 않는 격자선, 테두리, 배경을 제거합니다
- **의미 있게 정렬**: 자연적 순서(월, 단계)가 없다면 값 기준으로 카테고리를 정렬합니다
- **적절한 종횡비**: 시계열은 넓게(3:1에서 2:1), 비교는 더 정사각형에 가깝게
- **여백은 좋습니다**: 차트를 빽빽하게 채우지 마십시오. 각 시각화에 여유 공간을 줍니다

### 정확성

- **막대 차트는 0에서 시작**: 항상. 95에서 100까지의 막대는 5% 차이를 과장합니다
- **꺾은선 차트는 비제로 기준선 가능**: 변동 범위가 의미 있을 때
- **패널 간 일관된 스케일**: 여러 차트를 비교할 때 동일한 축 범위를 사용합니다
- **불확실성을 표시**: 데이터가 불확실할 때 오차 막대, 신뢰 구간 또는 범위를 사용합니다
- **축에 레이블을 부여**: 독자가 수치의 의미를 추측하게 하지 마십시오

## 접근성 고려사항

### 색맹

- 데이터 시리즈를 구분하는 데 색상에만 의존하지 않습니다
- 패턴 채우기, 다른 선 스타일(실선, 점선, 파선) 또는 직접 레이블을 추가합니다
- 색맹 시뮬레이터(예: Coblis, Sim Daltonism)로 테스트합니다
- 색맹 친화적 팔레트 사용: `sns.color_palette("colorblind")`

### 스크린 리더

- 차트의 핵심 발견을 설명하는 대체 텍스트를 포함합니다
- 시각화와 함께 데이터 테이블 대안을 제공합니다
- 시맨틱 제목과 레이블을 사용합니다

### 일반 접근성

- 데이터 요소와 배경 간 충분한 대비
- 레이블의 최소 텍스트 크기 10pt, 제목은 12pt
- 공간적 위치만으로 정보를 전달하는 것을 피합니다 (레이블을 추가)
- 인쇄를 고려합니다: 차트가 흑백에서도 작동하는가?

### 접근성 체크리스트

시각화를 공유하기 전에:
- [ ] 차트가 색상 없이도 작동 (패턴, 레이블 또는 선 스타일로 시리즈 구분)
- [ ] 표준 줌 수준에서 텍스트를 읽을 수 있음
- [ ] 제목이 데이터가 아닌 인사이트를 설명
- [ ] 축에 단위와 함께 레이블이 있음
- [ ] 범례가 명확하고 데이터를 가리지 않는 위치에 배치
- [ ] 데이터 소스와 날짜 범위가 기록됨
