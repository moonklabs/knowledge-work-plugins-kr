# Data Analyst 플러그인

Data Analyst 플러그인입니다. 주로 Anthropic의 agentic 데스크톱 애플리케이션인 [Cowork](https://claude.com/product/cowork)용으로 설계되었지만 Claude Code에서도 동작합니다. SQL query, data exploration, visualization, dashboard, insight generation을 지원합니다. 어떤 data warehouse, SQL dialect, analytics stack에서도 사용할 수 있습니다.

## 설치

```text
claude plugins add knowledge-work-plugins/data
```

## 주요 기능

이 플러그인은 Claude를 data analyst collaborator로 바꿉니다. Dataset exploration, optimized SQL 작성, visualization 생성, interactive dashboard 구축, stakeholder 공유 전 analysis validation을 돕습니다.

### Data Warehouse 연결이 있을 때

최상의 경험을 위해 data warehouse MCP server(Snowflake, Databricks, BigQuery 또는 SQL-compatible database)를 연결하세요. Claude는 다음을 수행합니다.

- Data warehouse를 직접 query합니다.
- Schema와 table metadata를 explore합니다.
- Copy/paste 없이 analysis를 end-to-end로 실행합니다.
- Result를 바탕으로 query를 iterate합니다.

### Data Warehouse 연결이 없을 때

Data warehouse 연결이 없으면 SQL result를 붙여 넣거나 CSV/Excel file을 upload해 analysis와 visualization을 수행합니다. Claude는 사용자가 수동 실행할 SQL query를 작성하고, 사용자가 제공한 result를 분석할 수도 있습니다.

## 명령

| 명령 | 설명 |
|---------|-------------|
| `/analyze` | Quick lookup부터 full analysis까지 data question에 답합니다. |
| `/explore-data` | Dataset의 shape, quality, pattern을 이해하기 위해 profile하고 explore합니다. |
| `/write-query` | 사용하는 dialect에 맞춰 best practice가 반영된 optimized SQL을 작성합니다. |
| `/create-viz` | Python으로 publication-quality visualization을 생성합니다. |
| `/build-dashboard` | Filter와 chart가 포함된 interactive HTML dashboard를 만듭니다. |
| `/validate` | 공유 전 methodology, accuracy, bias check로 analysis를 QA합니다. |

## 스킬

| 스킬 | 설명 |
|-------|-------------|
| `sql-queries` | 주요 data warehouse dialect(Snowflake, BigQuery, Databricks, PostgreSQL 등) 전반에서 정확하고 성능 좋은 SQL을 작성합니다. Query 작성, slow SQL optimization, dialect 간 translation, CTE/window function/aggregation이 있는 complex analytical query 작성에 사용합니다. |
| `data-exploration` | Data profiling, quality assessment, pattern discovery를 수행합니다. |
| `data-visualization` | Python(matplotlib, seaborn, plotly)으로 효과적인 data visualization을 만듭니다. Chart 작성, dataset에 맞는 chart type 선택, publication-quality figure 생성, accessibility와 color theory 같은 design principle 적용에 사용합니다. |
| `statistical-analysis` | Descriptive stats, trend analysis, outlier detection, hypothesis testing 등 statistical method를 적용합니다. Distribution analysis, significance testing, anomaly detection, correlation computation, statistical result interpretation에 사용합니다. |
| `data-validation` | Pre-delivery QA, sanity check, documentation standard를 다룹니다. |
| `interactive-dashboard-builder` | Chart.js, filter, styling을 사용해 HTML/JS dashboard를 구성합니다. |

## 예시 워크플로

### Ad-hoc analysis

```text
You: /analyze 지난 12개월 monthly revenue trend를 product line별로 보여줘.

Claude: [SQL query 작성] -> [Data warehouse에서 실행] -> [Trend chart 생성]
       -> [Key pattern 식별: "Product line A는 YoY 23% 성장했고 B는 flat"]
       -> [Sanity check로 result 검증]
```

### Data exploration

```text
You: /explore-data users table

Claude: [Table profile: 2.3M rows, 47 columns]
       -> [Report: created_at null 0.2%, email cardinality 99.8%]
       -> [Flag: status column에 unexpected value "UNKNOWN" 340 rows]
       -> [제안: "탐색 가치가 높은 dimension: plan_type, signup_source, country"]
```

### Query writing

```text
You: /write-query signup month별 cohort retention analysis가 필요해.
     1, 3, 6, 12개월 뒤 active 비율을 보여줘. Snowflake를 사용해.

Claude: [CTE가 포함된 optimized Snowflake SQL 작성]
       -> [각 step을 설명하는 comment 추가]
       -> [Partition pruning 관련 performance note 포함]
```

### Dashboard building

```text
You: /build-dashboard monthly revenue, top products,
     regional breakdown이 있는 sales dashboard를 만들어줘. Data는 여기 있어: [CSV 붙여넣기]

Claude: [Self-contained HTML file 생성]
       -> [Interactive Chart.js visualization 포함]
       -> [Region 및 time period dropdown filter 추가]
       -> [Review를 위해 browser에서 열기]
```

### 공유 전 validation

```text
You: /validate [shares analysis document]

Claude: [Methodology review] -> [Churn analysis의 survivorship bias 확인]
       -> [Aggregation logic 검증] -> [Flag: "Denominator가 trial user를 제외해
          conversion rate를 약 5pp 과대평가할 수 있음"]
       -> [Confidence: "Caveat를 명시하면 공유 가능"]
```

## Data Stack 연결

> 익숙하지 않은 placeholder가 보이거나 연결된 도구를 확인해야 한다면 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

이 플러그인은 data infrastructure에 연결될 때 가장 잘 동작합니다. 다음 MCP server를 추가하세요.

- **Data Warehouse**: Snowflake, Databricks, BigQuery, Definite 또는 SQL-compatible database
- **Analytics/BI**: Amplitude, Looker, Tableau 등
- **Notebooks**: Jupyter, Hex 등
- **Spreadsheets**: Google Sheets, Excel
- **Data Orchestration**: Airflow, dbt, Dagster, Prefect
- **Data Ingestion**: Fivetran, Airbyte, Stitch

Direct data access를 활성화하려면 `.mcp.json` 또는 Claude Code settings에서 MCP server를 설정하세요.
