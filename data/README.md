# Data Analyst 플러그인

[Cowork](https://claude.com/product/cowork)를 위한 데이터 분석가 플러그인입니다. Anthropic의 에이전트형 데스크톱 앱인 Cowork에 맞춰 설계했지만 Claude Code에서도 동작합니다. SQL 쿼리, 데이터 탐색, 시각화, 대시보드, 인사이트 생성까지 지원합니다. 어떤 데이터 웨어하우스, 어떤 SQL 방언, 어떤 분석 스택에서도 동작합니다.

## 설치

```
claude plugins add knowledge-work-plugins/data
```

## 제공 기능

이 플러그인은 Claude를 데이터 분석 협업자로 바꿉니다. 데이터셋 탐색, 최적화된 SQL 작성, 시각화 생성, 인터랙티브 대시보드 제작, 이해관계자 공유 전 분석 검증을 도와줍니다.

### 데이터 웨어하우스 연결 시

최상의 경험을 위해 데이터 웨어하우스 MCP 서버(예: Snowflake, Databricks, BigQuery 또는 SQL 호환 데이터베이스)를 연결하세요. Claude는 다음을 수행합니다.

- 데이터 웨어하우스를 직접 쿼리
- 스키마 및 테이블 메타데이터 탐색
- 복사/붙여넣기 없이 분석을 엔드투엔드로 실행
- 결과에 맞춰 쿼리를 반복 개선

### 데이터 웨어하우스 미연결 시

데이터 웨어하우스가 없으면 SQL 결과를 붙여넣거나 CSV/Excel 파일을 업로드해 분석과 시각화를 수행합니다. Claude가 SQL 쿼리를 작성해 주면 사용자가 직접 실행한 뒤 결과를 제공해 분석할 수도 있습니다.

## 커맨드

| Command | 설명 |
|---------|-------------|
| `/analyze` | 데이터 질문에 답변: 빠른 조회부터 전체 분석까지 |
| `/explore-data` | 데이터셋의 형태, 품질, 패턴을 이해하기 위한 프로파일링/탐색 |
| `/write-query` | 방언별 모범 사례를 반영한 최적화 SQL 작성 |
| `/create-viz` | Python으로 출판 수준의 시각화 생성 |
| `/build-dashboard` | 필터와 차트가 있는 인터랙티브 HTML 대시보드 생성 |
| `/validate` | 공유 전 QA: 방법론, 정확성, 편향 점검 |

## 스킬

| Skill | 설명 |
|-------|-------------|
| `sql-queries` | 방언별 SQL 모범 사례, 공통 패턴, 성능 최적화 |
| `data-exploration` | 데이터 프로파일링, 품질 평가, 패턴 탐색 |
| `data-visualization` | 차트 선택, Python 시각화 코드 패턴, 디자인 원칙 |
| `statistical-analysis` | 기술 통계, 추세 분석, 이상치 탐지, 가설 검정 |
| `data-validation` | 전달 전 QA, 합리성 체크, 문서화 기준 |
| `interactive-dashboard-builder` | Chart.js 기반 HTML/JS 대시보드 구성, 필터, 스타일링 |

## 예시 워크플로

### 애드혹 분석

```
You: /analyze 지난 12개월 동안 제품 라인별 월간 매출 추세는 어땠나요?

Claude: [SQL 쿼리 작성] → [데이터 웨어하우스에서 실행] → [추세 차트 생성]
       → [핵심 패턴 식별: "제품 라인 A는 YoY 23% 성장, B는 정체"]
       → [합리성 체크로 결과 검증]
```

### 데이터 탐색

```
You: /explore-data users 테이블

Claude: [테이블 프로파일: 2.3M 행, 47개 컬럼]
       → [보고: created_at 결측 0.2%, email 카디널리티 99.8%]
       → [플래그: status 컬럼에 예상치 않은 값 "UNKNOWN" 340건]
       → [제안: "우선 탐색할 고가치 차원: plan_type, signup_source, country"]
```

### 쿼리 작성

```
You: /write-query 가입 월 기준 코호트 잔존 분석이 필요합니다.
     1, 3, 6, 12개월 후에도 활성 사용자 비율을 보여 주세요. Snowflake를 사용합니다.

Claude: [CTE를 사용한 최적화 Snowflake SQL 작성]
       → [각 단계 설명 주석 추가]
       → [파티션 프루닝에 대한 성능 메모 포함]
```

### 대시보드 구축

```
You: /build-dashboard 월간 매출, 상위 제품, 지역별 분해를 포함한
     세일즈 대시보드를 만들어 주세요. 데이터는 여기 있어요: [CSV 붙여넣기]

Claude: [자체 포함형 HTML 파일 생성]
       → [Chart.js 기반 인터랙티브 시각화 포함]
       → [지역/기간 드롭다운 필터 추가]
       → [브라우저에서 열어 리뷰]
```

### 공유 전 검증

```
You: /validate [분석 문서 공유]

Claude: [방법론 리뷰] → [이탈 분석에서 생존자 편향 확인]
       → [집계 로직 검증] → [플래그: "분모에서 체험판 사용자가 제외되어
          전환율이 약 5pp 과대 추정될 수 있음"]
       → [신뢰도: "주의사항을 포함해 공유 가능"]
```

## 데이터 스택 연결

> 낯선 플레이스홀더가 보이거나 어떤 도구가 연결되어 있는지 확인하려면 [CONNECTORS.md](CONNECTORS.md)를 보세요.

이 플러그인은 데이터 인프라에 연결할 때 가장 효과적입니다. 다음 MCP 서버를 추가하세요.

- **Data Warehouse**: Snowflake, Databricks, BigQuery 또는 SQL 호환 데이터베이스
- **Analytics/BI**: Amplitude, Looker, Tableau 등
- **Notebooks**: Jupyter, Hex 등
- **Spreadsheets**: Google Sheets, Excel
- **Data Orchestration**: Airflow, dbt, Dagster, Prefect
- **Data Ingestion**: Fivetran, Airbyte, Stitch

`.mcp.json` 또는 Claude Code 설정에서 MCP 서버를 구성해 직접 데이터 접근을 활성화하세요.
