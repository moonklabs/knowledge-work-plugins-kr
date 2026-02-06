# Finance & Accounting 플러그인

[Cowork](https://claude.com/product/cowork)를 위한 재무/회계 플러그인입니다. Anthropic의 에이전트형 데스크톱 앱인 Cowork에 맞춰 설계했지만 Claude Code에서도 동작합니다. 월말 마감, 전표 준비, 계정 조정, 재무제표 생성, 변동 분석, SOX 감사 지원을 제공합니다.

> **중요:** 이 플러그인은 재무/회계 워크플로를 보조하지만 재무, 세무, 감사 자문을 제공하지 않습니다. 모든 결과물은 재무 보고, 규제 제출, 감사 문서에 사용하기 전에 자격을 갖춘 전문가의 검토가 필요합니다.

## 설치

```bash
claude plugins add knowledge-work-plugins/finance
```

## 커맨드

| Command | 설명 |
|---------|-------------|
| `/journal-entry` | 전표 준비: 적절한 차/대변과 증빙을 포함해 발생주의, 고정자산, 선급비용, 급여, 수익 전표 생성 |
| `/reconciliation` | 계정 조정: GL 잔액을 서브레저, 은행, 또는 제3자 잔액과 비교하고 조정 항목 식별 |
| `/income-statement` | 손익계산서 생성: 전기 대비 비교 및 변동 분석 포함 |
| `/variance-analysis` | 변동/플럭스 분석: 변동 요인을 분해하고 서술형 설명과 워터폴 분석 제공 |
| `/sox-testing` | SOX 컴플라이언스 테스트: 샘플 선정, 테스트 워크페이퍼, 통제 평가 생성 |

## 스킬

| Skill | 설명 |
|-------|-------------|
| `journal-entry-prep` | JE 준비 모범 사례, 표준 발생주의 유형, 증빙 요건, 리뷰 워크플로 |
| `reconciliation` | GL-서브레저, 은행 조정, 계열사 간 조정 방법론과 조정 항목 분류/노후화 |
| `financial-statements` | 손익계산서, 대차대조표, 현금흐름표 포맷과 GAAP 표시/플럭스 분석 방법론 |
| `variance-analysis` | 변동 분해 기법(가격/수량, 요율/믹스), 중요성 기준, 서술 생성, 워터폴 차트 |
| `close-management` | 월말 마감 체크리스트, 작업 순서, 의존성, 상태 추적, 일자별 활동 |
| `audit-support` | SOX 404 통제 테스트 방법론, 샘플 선정, 문서화 기준, 결함 분류 |

## 예시 워크플로

### 월말 마감

1. `/journal-entry ap-accrual 2024-12` 실행으로 매입채무 발생주의 전표 생성
2. `/journal-entry prepaid 2024-12` 실행으로 선급비용 상각
3. `/journal-entry fixed-assets 2024-12` 실행으로 감가상각 계상
4. `/reconciliation cash 2024-12` 실행으로 은행 계정 조정
5. `/reconciliation accounts-receivable 2024-12` 실행으로 매출채권 서브레저 조정
6. `/income-statement monthly 2024-12` 실행으로 플럭스 분석 포함 P&L 생성

### 변동 분석

1. `/variance-analysis revenue 2024-Q4 vs 2024-Q3` 실행으로 매출 변동 분석
2. `/variance-analysis opex 2024-12 vs budget` 실행으로 판관비 변동 조사
3. 워터폴 분석을 검토하고 설명되지 않은 변동에 대한 컨텍스트 제공

### SOX 테스트

1. `/sox-testing revenue-recognition 2024-Q4` 실행으로 매출 인식 통제 테스트 워크페이퍼 생성
2. `/sox-testing procure-to-pay 2024-Q4` 실행으로 구매 통제 테스트
3. 샘플 선정 결과를 검토하고 테스트 결과 문서화

## MCP 통합

> 낯선 플레이스홀더가 보이거나 어떤 도구가 연결되어 있는지 확인하려면 [CONNECTORS.md](CONNECTORS.md)를 보세요.

이 플러그인은 MCP 서버로 재무 데이터 소스에 연결할 때 가장 효과적입니다. 관련 서버를 `.mcp.json`에 추가하세요.

### ERP / 회계 시스템

ERP(예: NetSuite, SAP) MCP 서버를 연결해 시산표, 서브레저 데이터, 전표를 자동으로 가져옵니다.

### 데이터 웨어하우스

데이터 웨어하우스(예: Snowflake, BigQuery) MCP 서버를 연결해 재무 데이터를 조회하고 변동 분석과 과거 비교를 수행합니다.

### 스프레드시트

스프레드시트 도구(예: Google Sheets, Excel)를 연결해 워크페이퍼 생성, 조정 템플릿, 재무 모델 업데이트를 수행합니다.

### 분석 / BI

BI 플랫폼(예: Tableau, Looker)을 연결해 대시보드, KPI, 추세 데이터를 가져와 변동 설명에 활용합니다.

> **참고:** ERP와 데이터 웨어하우스 MCP 서버를 연결해 재무 데이터를 자동으로 가져오세요. 연결하지 않으면 데이터를 붙여넣거나 파일을 업로드해 분석할 수 있습니다.

## 설정

이 플러그인 디렉터리의 `.mcp.json`에 `mcpServers` 섹션으로 데이터 소스 MCP 서버를 추가하세요. `recommendedCategories` 필드는 플러그인 기능을 강화하는 통합 유형을 나열합니다.

- `erp-accounting` — GL, 서브레저, JE 데이터를 위한 ERP 또는 회계 시스템
- `data-warehouse` — 재무 쿼리와 과거 데이터 분석을 위한 데이터 웨어하우스
- `spreadsheets` — 워크페이퍼 생성을 위한 스프레드시트 도구
- `analytics-bi` — 대시보드와 KPI 데이터를 위한 BI 도구
- `documents` — 정책, 메모, 증빙 문서를 위한 문서 저장소
- `email` — 보고서 발송 및 승인 요청을 위한 이메일
- `chat` — 마감 진행 상황 공유 및 질문을 위한 팀 커뮤니케이션
