# Finance & Accounting 플러그인

Finance 및 accounting plugin입니다. 주로 Anthropic의 agentic desktop application인 [Cowork](https://claude.com/product/cowork)용으로 설계되었지만 Claude Code에서도 동작합니다. Month-end close, journal entry preparation, account reconciliation, financial statement generation, variance analysis, SOX audit support를 지원합니다.

> **중요:** 이 플러그인은 finance 및 accounting workflow를 보조하지만 financial, tax, audit advice를 제공하지 않습니다. 모든 output은 financial reporting, regulatory filing, audit documentation에 사용하기 전에 자격을 갖춘 financial professional의 검토를 받아야 합니다.

## 설치

```bash
claude plugins add knowledge-work-plugins/finance
```

## 명령

| 명령 | 설명 |
|---------|-------------|
| `/journal-entry` | 적절한 debit/credit과 supporting detail이 포함된 accrual, fixed asset, prepaid, payroll, revenue journal entry를 준비합니다. |
| `/reconciliation` | GL balance를 subledger, bank, third-party balance와 비교하고 reconciling item을 식별합니다. |
| `/income-statement` | Period-over-period comparison과 variance analysis가 포함된 P&L을 생성합니다. |
| `/variance-analysis` | Variance를 driver로 분해하고 narrative explanation과 waterfall analysis를 제공합니다. |
| `/sox-testing` | Sample selection, testing workpaper, control assessment를 생성해 SOX compliance testing을 지원합니다. |

## 스킬

| 스킬 | 설명 |
|-------|-------------|
| `journal-entry-prep` | Month-end close를 위해 적절한 debit, credit, supporting documentation을 갖춘 journal entry를 준비합니다. Accrual, prepaid amortization, fixed asset depreciation, payroll entry, revenue recognition 또는 manual journal entry booking에 사용합니다. |
| `reconciliation` | GL balance를 subledger, bank statement, third-party data와 비교해 account를 reconcile합니다. Bank reconciliation, GL-to-subledger rec, intercompany reconciliation, reconciling item 식별 및 categorization에 사용합니다. |
| `financial-statements` | Period-over-period comparison과 variance analysis가 포함된 financial statement(income statement, balance sheet, cash flow)를 생성합니다. Monthly/quarterly P&L 준비, book closing 중 material variance flagging, actual vs budget comparison, leadership review용 financial summary 작성, GAAP presentation requirement와 period-end adjustment 확인에 사용합니다. |
| `variance-analysis` | Financial variance를 narrative explanation과 waterfall analysis가 포함된 driver로 분해합니다. Budget vs actual, period-over-period change, revenue/expense variance 분석, leadership용 variance commentary 준비에 사용합니다. |
| `close-management` | Task sequencing, dependency, status tracking으로 month-end close process를 관리합니다. Close calendar planning, close progress tracking, blocker identification, close activity day-by-day sequencing에 사용합니다. |
| `audit-support` | Control testing methodology, sample selection, documentation standard로 SOX 404 compliance를 지원합니다. Testing workpaper 생성, audit sample 선택, control deficiency classification, internal/external audit 준비에 사용합니다. |

## 예시 워크플로

### Month-end close

1. `/journal-entry ap-accrual 2024-12`를 실행해 AP accrual entry를 생성합니다.
2. `/journal-entry prepaid 2024-12`를 실행해 prepaid expense를 amortize합니다.
3. `/journal-entry fixed-assets 2024-12`를 실행해 depreciation을 book합니다.
4. `/reconciliation cash 2024-12`를 실행해 bank account를 reconcile합니다.
5. `/reconciliation accounts-receivable 2024-12`를 실행해 AR subledger를 reconcile합니다.
6. `/income-statement monthly 2024-12`를 실행해 flux analysis가 포함된 P&L을 생성합니다.

### Variance analysis

1. `/variance-analysis revenue 2024-Q4 vs 2024-Q3`를 실행해 revenue variance를 분석합니다.
2. `/variance-analysis opex 2024-12 vs budget`을 실행해 operating expense variance를 조사합니다.
3. Waterfall analysis를 review하고 unexplained variance에 대한 context를 제공합니다.

### SOX testing

1. `/sox-testing revenue-recognition 2024-Q4`를 실행해 revenue control testing workpaper를 생성합니다.
2. `/sox-testing procure-to-pay 2024-Q4`를 실행해 procurement control을 test합니다.
3. Sample selection을 review하고 test result를 문서화합니다.

## MCP 통합

> 익숙하지 않은 placeholder가 보이거나 연결된 tool을 확인해야 한다면 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

이 plugin은 MCP server를 통해 financial data source에 연결할 때 가장 잘 동작합니다. 관련 server를 `.mcp.json`에 추가하세요.

### ERP / accounting system

ERP(예: NetSuite, SAP) MCP server를 연결해 trial balance, subledger data, journal entry를 자동으로 가져옵니다.

### Data warehouse

Data warehouse(예: Snowflake, BigQuery) MCP server를 연결해 financial data를 query하고 variance analysis를 실행하며 historical comparison을 가져옵니다.

### Spreadsheet

Workpaper generation, reconciliation template, financial model update를 위해 spreadsheet tool(예: Google Sheets, Excel)을 연결합니다.

### Analytics / BI

BI platform(예: Tableau, Looker)을 연결해 variance explanation에 필요한 dashboard, KPI, trend data를 가져옵니다.

> **Note:** Financial data를 자동으로 가져오려면 ERP 및 data warehouse MCP server를 연결하세요. 없으면 data를 paste하거나 file을 upload해 analysis할 수 있습니다.

## 설정

이 plugin directory의 `.mcp.json`에서 `mcpServers` section에 data source MCP server를 추가하세요. `recommendedCategories` field는 plugin 기능을 강화하는 integration type을 나열합니다.

- `erp-accounting` — GL, subledger, JE data용 ERP 또는 accounting system
- `data-warehouse` — Financial query 및 historical data용 data warehouse
- `spreadsheets` — Workpaper generation용 spreadsheet tool
- `analytics-bi` — Dashboard 및 KPI data용 BI tool
- `documents` — Policy, memo, support 자료용 document storage
- `email` — Report 발송 및 approval request용 email
- `chat` — Close status update 및 question용 team communication
