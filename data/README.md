# Data Analyst 플러그인

Data analyst plugin입니다. 주로 Anthropic의 agentic desktop application인 [Cowork](https://claude.com/product/cowork)용으로 설계되었지만 Claude Code에서도 동작합니다. SQL query, data exploration, visualization, dashboard, insight generation을 지원합니다. 어떤 data warehouse, SQL dialect, analytics stack에서도 사용할 수 있습니다.

## 설치

```
claude plugins add knowledge-work-plugins/data
```

## 주요 기능

This plugin transforms Claude into a data analyst collaborator. It helps you explore datasets, write optimized SQL, build visualizations, create interactive dashboards, and validate analyses before sharing with stakeholders.

### With a Data Warehouse Connection

Connect your data warehouse MCP server (e.g., Snowflake, Databricks, BigQuery, or any SQL-compatible database) for the best experience. Claude will:

- Query your data warehouse directly
- Explore schemas and table metadata
- Run analyses end-to-end without copy-pasting
- Iterate on queries based on results

### Without a Data Warehouse Connection

Without a data warehouse connection, paste SQL results or upload CSV/Excel files for analysis and visualization. Claude can also write SQL queries for you to run manually, and then analyze the results you provide.

## 명령

| 명령 | 설명 |
|---------|-------------|
| `/analyze` | Answer data questions -- from quick lookups to full analyses |
| `/explore-data` | Profile and explore a dataset to understand its shape, quality, and patterns |
| `/write-query` | Write optimized SQL for your dialect with best practices |
| `/create-viz` | Create publication-quality visualizations with Python |
| `/build-dashboard` | Build interactive HTML dashboards with filters and charts |
| `/validate` | QA an analysis before sharing -- methodology, accuracy, and bias checks |

## 스킬

| 스킬 | 설명 |
|-------|-------------|
| `sql-queries` | 주요 data warehouse dialect(Snowflake, BigQuery, Databricks, PostgreSQL 등) 전반에서 correct하고 performant한 SQL을 작성합니다. Query 작성, slow SQL optimization, dialect 간 translation, CTE/window function/aggregation이 있는 complex analytical query 작성에 사용합니다. |
| `data-exploration` | Data profiling, quality assessment, and pattern discovery |
| `data-visualization` | Python(matplotlib, seaborn, plotly)으로 효과적인 data visualization을 만듭니다. Chart 작성, dataset에 맞는 chart type 선택, publication-quality figure 생성, accessibility와 color theory 같은 design principle 적용에 사용합니다. |
| `statistical-analysis` | Descriptive stats, trend analysis, outlier detection, hypothesis testing 등 statistical method를 적용합니다. Distribution analysis, significance testing, anomaly detection, correlation computation, statistical result interpretation에 사용합니다. |
| `data-validation` | Pre-delivery QA, sanity checks, and documentation standards |
| `interactive-dashboard-builder` | HTML/JS dashboard construction with Chart.js, filters, and styling |

## 예시 워크플로

### Ad-Hoc Analysis

```
You: /analyze What was our monthly revenue trend for the past 12 months, broken down by product line?

Claude: [Writes SQL query] → [Executes against data warehouse] → [Generates trend chart]
       → [Identifies key patterns: "Product line A grew 23% YoY while B was flat"]
       → [Validates results with sanity checks]
```

### Data Exploration

```
You: /explore-data users table

Claude: [Profiles table: 2.3M rows, 47 columns]
       → [Reports: created_at has 0.2% nulls, email has 99.8% cardinality]
       → [Flags: status column has unexpected value "UNKNOWN" in 340 rows]
       → [Suggests: "High-value dimensions to explore: plan_type, signup_source, country"]
```

### Query Writing

```
You: /write-query I need a cohort retention analysis -- users grouped by signup month,
     showing what % are still active 1, 3, 6, and 12 months later. We use Snowflake.

Claude: [Writes optimized Snowflake SQL with CTEs]
       → [Adds comments explaining each step]
       → [Includes performance notes about partition pruning]
```

### Dashboard Building

```
You: /build-dashboard Create a sales dashboard with monthly revenue, top products,
     and regional breakdown. Here's the data: [pastes CSV]

Claude: [Generates self-contained HTML file]
       → [Includes interactive Chart.js visualizations]
       → [Adds dropdown filters for region and time period]
       → [Opens in browser for review]
```

### Pre-Share Validation

```
You: /validate [shares analysis document]

Claude: [Reviews methodology] → [Checks for survivorship bias in churn analysis]
       → [Verifies aggregation logic] → [Flags: "Denominator excludes trial users
          which could overstate conversion rate by ~5pp"]
       → [Confidence: "Ready to share with noted caveat"]
```

## Data Stack 연결

> If you see unfamiliar placeholders or need to check which tools are connected, see [CONNECTORS.md](CONNECTORS.md).

This plugin works best when connected to your data infrastructure. Add MCP servers for:

- **Data Warehouse**: Snowflake, Databricks, BigQuery, Definite, or any SQL-compatible database
- **Analytics/BI**: Amplitude, Looker, Tableau, or similar
- **Notebooks**: Jupyter, Hex, or similar
- **Spreadsheets**: Google Sheets, Excel
- **Data Orchestration**: Airflow, dbt, Dagster, Prefect
- **Data Ingestion**: Fivetran, Airbyte, Stitch

Configure MCP servers in your `.mcp.json` or Claude Code settings to enable direct data access.
