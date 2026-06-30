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
| `/journal-entry` | Journal entry preparation — generate accruals, fixed asset entries, prepaids, payroll, and revenue entries with proper debits/credits and supporting detail |
| `/reconciliation` | Account reconciliation — compare GL balances to subledger, bank, or third-party balances and identify reconciling items |
| `/income-statement` | Income statement generation — produce P&L with period-over-period comparison and variance analysis |
| `/variance-analysis` | Variance/flux analysis — decompose variances into drivers with narrative explanations and waterfall analysis |
| `/sox-testing` | SOX compliance testing — generate sample selections, testing workpapers, and control assessments |

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

### Month-End Close

1. Run `/journal-entry ap-accrual 2024-12` to generate AP accrual entries
2. Run `/journal-entry prepaid 2024-12` to amortize prepaid expenses
3. Run `/journal-entry fixed-assets 2024-12` to book depreciation
4. Run `/reconciliation cash 2024-12` to reconcile bank accounts
5. Run `/reconciliation accounts-receivable 2024-12` to reconcile AR subledger
6. Run `/income-statement monthly 2024-12` to generate the P&L with flux analysis

### Variance Analysis

1. Run `/variance-analysis revenue 2024-Q4 vs 2024-Q3` to analyze revenue variances
2. Run `/variance-analysis opex 2024-12 vs budget` to investigate operating expense variances
3. Review the waterfall analysis and provide context on unexplained variances

### SOX Testing

1. Run `/sox-testing revenue-recognition 2024-Q4` to generate revenue control testing workpapers
2. Run `/sox-testing procure-to-pay 2024-Q4` to test procurement controls
3. Review sample selections and document test results

## MCP Integration

> If you see unfamiliar placeholders or need to check which tools are connected, see [CONNECTORS.md](CONNECTORS.md).

This plugin works best when connected to your financial data sources via MCP servers. Add the relevant servers to your `.mcp.json`:

### ERP / Accounting System

Connect your ERP (e.g., NetSuite, SAP) MCP server to pull trial balances, subledger data, and journal entries automatically.

### Data Warehouse

Connect your data warehouse (e.g., Snowflake, BigQuery) MCP server to query financial data, run variance analysis, and pull historical comparisons.

### Spreadsheets

Connect spreadsheet tools (e.g., Google Sheets, Excel) for workpaper generation, reconciliation templates, and financial model updates.

### Analytics / BI

Connect your BI platform (e.g., Tableau, Looker) to pull dashboards, KPIs, and trend data for variance explanations.

> **Note:** Connect your ERP and data warehouse MCP servers to pull financial data automatically. Without these, you can paste data or upload files for analysis.

## Configuration

Add your data source MCP servers to the `mcpServers` section of `.mcp.json` in this plugin directory. The `recommendedCategories` field lists the types of integrations that enhance this plugin's capabilities:

- `erp-accounting` — ERP or accounting system for GL, subledger, and JE data
- `data-warehouse` — Data warehouse for financial queries and historical data
- `spreadsheets` — Spreadsheet tools for workpaper generation
- `analytics-bi` — BI tools for dashboards and KPI data
- `documents` — Document storage for policies, memos, and support
- `email` — Email for sending reports and requesting approvals
- `chat` — Team communication for close status updates and questions
