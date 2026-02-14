# 재무 및 회계 플러그인

Anthropic의 에이전트형 데스크톱 애플리케이션인 [Cowork](https://claude.com/product/cowork)을 위해 설계된 재무 및 회계 플러그인입니다. Claude Code에서도 사용할 수 있습니다. 월마감(month-end close), 분개(journal entry) 준비, 계정 조정(reconciliation), 재무제표 생성, 차이 분석(variance analysis), SOX 감사 지원을 제공합니다.

> **중요**: 이 플러그인은 재무 및 회계 워크플로우를 지원하지만 재무, 세무 또는 감사 자문을 제공하지는 않습니다. 모든 결과물은 재무 보고, 규제 제출 또는 감사 문서에 사용하기 전에 자격을 갖춘 재무 전문가의 검토를 받아야 합니다.

## 설치

```bash
claude plugins add knowledge-work-plugins/finance
```

## 커맨드

| 커맨드 | 설명 |
|---------|------|
| `/journal-entry` | 분개 준비 — 발생, 고정자산, 선급비용, 급여, 수익 항목에 대해 적절한 차변/대변 및 증빙 자료와 함께 분개를 생성합니다 |
| `/reconciliation` | 계정 조정 — 총계정원장 잔액을 보조원장, 은행 또는 제3자 잔액과 비교하고 조정 항목을 식별합니다 |
| `/income-statement` | 손익계산서 생성 — 기간별 비교 및 차이 분석이 포함된 손익계산서(P&L)를 작성합니다 |
| `/variance-analysis` | 차이/변동 분석 — 차이를 동인별로 분해하고 서술적 설명 및 워터폴 분석을 제공합니다 |
| `/sox-testing` | SOX 컴플라이언스 테스트 — 표본 선정, 테스트 조서 및 통제 평가를 생성합니다 |

## 스킬

| 스킬 | 설명 |
|------|------|
| `journal-entry-prep` | 분개 준비 모범 사례, 표준 발생 유형, 증빙 문서 요건 및 검토 워크플로우 |
| `reconciliation` | 총계정원장-보조원장, 은행 조정, 회사간 조정을 위한 조정 방법론 및 조정 항목 분류와 경과 기간 분석 |
| `financial-statements` | GAAP 표시 기준과 변동 분석 방법론이 포함된 손익계산서, 재무상태표, 현금흐름표 형식 |
| `variance-analysis` | 차이 분해 기법(가격/수량, 비율/믹스), 중요성 기준, 서술 생성 및 워터폴 차트 |
| `close-management` | 월마감 체크리스트, 작업 순서, 의존관계, 진행 상태 추적 및 일별 주요 마감 활동 |
| `audit-support` | SOX 404 통제 테스트 방법론, 표본 선정, 문서화 기준 및 결함 분류 |

## 예시 워크플로우

### 월마감

1. `/journal-entry ap-accrual 2024-12`를 실행하여 매입채무 발생 분개를 생성합니다
2. `/journal-entry prepaid 2024-12`를 실행하여 선급비용을 상각합니다
3. `/journal-entry fixed-assets 2024-12`를 실행하여 감가상각을 기장합니다
4. `/reconciliation cash 2024-12`를 실행하여 은행 계좌를 조정합니다
5. `/reconciliation accounts-receivable 2024-12`를 실행하여 매출채권 보조원장을 조정합니다
6. `/income-statement monthly 2024-12`를 실행하여 변동 분석이 포함된 손익계산서를 생성합니다

### 차이 분석

1. `/variance-analysis revenue 2024-Q4 vs 2024-Q3`를 실행하여 수익 차이를 분석합니다
2. `/variance-analysis opex 2024-12 vs budget`를 실행하여 영업비용 차이를 조사합니다
3. 워터폴 분석을 검토하고 설명되지 않은 차이에 대한 맥락을 제공합니다

### SOX 테스트

1. `/sox-testing revenue-recognition 2024-Q4`를 실행하여 수익 인식 통제 테스트 조서를 생성합니다
2. `/sox-testing procure-to-pay 2024-Q4`를 실행하여 구매 통제를 테스트합니다
3. 표본 선정을 검토하고 테스트 결과를 문서화합니다

## MCP 연동

> 익숙하지 않은 플레이스홀더가 보이거나 연결된 도구를 확인해야 하는 경우, [CONNECTORS.md](CONNECTORS.md)를 참조하십시오.

이 플러그인은 MCP 서버를 통해 재무 데이터 소스에 연결할 때 가장 효과적으로 동작합니다. `.mcp.json`에 관련 서버를 추가하십시오.

### ERP / 회계 시스템

ERP(예: NetSuite, SAP) MCP 서버를 연결하면 시산표, 보조원장 데이터, 분개를 자동으로 가져올 수 있습니다.

### 데이터 웨어하우스

데이터 웨어하우스(예: Snowflake, BigQuery) MCP 서버를 연결하면 재무 데이터 조회, 차이 분석 실행, 과거 비교 데이터 추출이 가능합니다.

### 스프레드시트

스프레드시트 도구(예: Google Sheets, Excel)를 연결하면 조서 생성, 조정 템플릿, 재무 모델 업데이트에 활용할 수 있습니다.

### 분석 / BI

BI 플랫폼(예: Tableau, Looker)을 연결하면 대시보드, KPI, 추세 데이터를 차이 설명에 활용할 수 있습니다.

> **참고:** ERP 및 데이터 웨어하우스 MCP 서버를 연결하면 재무 데이터를 자동으로 가져올 수 있습니다. 연결하지 않아도 데이터를 붙여넣기하거나 파일을 업로드하여 분석할 수 있습니다.

## 구성

이 플러그인 디렉토리의 `.mcp.json` 내 `mcpServers` 섹션에 데이터 소스 MCP 서버를 추가하십시오. `recommendedCategories` 필드에는 이 플러그인의 기능을 향상시키는 연동 유형이 나열되어 있습니다.

- `erp-accounting` — 총계정원장, 보조원장, 분개 데이터를 위한 ERP 또는 회계 시스템
- `data-warehouse` — 재무 조회 및 과거 데이터를 위한 데이터 웨어하우스
- `spreadsheets` — 조서 생성을 위한 스프레드시트 도구
- `analytics-bi` — 대시보드 및 KPI 데이터를 위한 BI 도구
- `documents` — 정책, 메모, 증빙 자료를 위한 문서 저장소
- `email` — 보고서 발송 및 승인 요청을 위한 이메일
- `chat` — 마감 진행 상황 업데이트 및 질의를 위한 팀 커뮤니케이션
