# Data Analyst 플러그인

Anthropic의 에이전트 데스크톱 애플리케이션인 [Cowork](https://claude.com/product/cowork)을 위해 설계된 데이터 분석 플러그인입니다. Claude Code에서도 사용할 수 있습니다. SQL 쿼리, 데이터 탐색, 시각화, 대시보드, 인사이트 생성 기능을 제공합니다. 모든 데이터 웨어하우스, SQL 방언, 분석 스택에서 작동합니다.

## 설치

```
claude plugins add knowledge-work-plugins/data
```

## 기능 소개

이 플러그인은 Claude를 데이터 분석 협업 파트너로 변환합니다. 데이터셋 탐색, 최적화된 SQL 작성, 시각화 구축, 대화형 대시보드 생성, 이해관계자에게 공유하기 전 분석 검증을 지원합니다.

### 데이터 웨어하우스 연결 시

데이터 웨어하우스 MCP 서버(예: Snowflake, Databricks, BigQuery 또는 기타 SQL 호환 데이터베이스)를 연결하면 최상의 경험을 제공합니다. Claude가 다음을 수행합니다:

- 데이터 웨어하우스에 직접 쿼리 실행
- 스키마 및 테이블 메타데이터 탐색
- 복사/붙여넣기 없이 분석을 처음부터 끝까지 실행
- 결과를 기반으로 쿼리 반복 개선

### 데이터 웨어하우스 미연결 시

데이터 웨어하우스 연결 없이도 SQL 결과를 붙여넣거나 CSV/Excel 파일을 업로드하여 분석 및 시각화가 가능합니다. Claude가 수동으로 실행할 SQL 쿼리를 작성해 주고, 제공된 결과를 분석할 수도 있습니다.

## 커맨드

| 커맨드 | 설명 |
|---------|------|
| `/analyze` | 데이터 질문에 답변 -- 간단한 조회부터 전체 분석까지 |
| `/explore-data` | 데이터셋의 형태, 품질, 패턴을 프로파일링하고 탐색 |
| `/write-query` | 모범 사례를 적용한 SQL 방언별 최적화 쿼리 작성 |
| `/create-viz` | Python으로 출판 수준의 시각화 생성 |
| `/build-dashboard` | 필터와 차트가 포함된 대화형 HTML 대시보드 구축 |
| `/validate` | 공유 전 분석 QA -- 방법론, 정확성, 편향 검사 |

## 스킬

| 스킬 | 설명 |
|------|------|
| `sql-queries` | 방언별 SQL 모범 사례, 일반적인 패턴, 성능 최적화 |
| `data-exploration` | 데이터 프로파일링, 품질 평가, 패턴 발견 |
| `data-visualization` | 차트 선택, Python 시각화 코드 패턴, 디자인 원칙 |
| `statistical-analysis` | 기술 통계, 추세 분석, 이상치 탐지, 가설 검정 |
| `data-validation` | 전달 전 QA, 정합성 검사, 문서화 표준 |
| `interactive-dashboard-builder` | Chart.js, 필터, 스타일링을 활용한 HTML/JS 대시보드 구축 |

## 워크플로우 예시

### 애드혹 분석

```
나: /analyze 지난 12개월 동안 제품 라인별 월간 수익 추세는 어땠나요?

Claude: [SQL 쿼리 작성] → [데이터 웨어하우스에 대해 실행] → [추세 차트 생성]
       → [핵심 패턴 식별: "제품 라인 A는 전년 대비 23% 성장한 반면 B는 정체됨"]
       → [합리성 검사로 결과 검증]
```

### 데이터 탐색

```
나: /explore-data users 테이블

Claude: [테이블 프로파일링: 230만 행, 47개 열]
       → [보고: created_at에 0.2%의 null이 있음, email은 99.8%의 카디널리티를 가짐]
       → [플래그: status 열에 340개 행에서 예상치 못한 값 "UNKNOWN"이 있음]
       → [제안: "탐색할 고가치 차원: plan_type, signup_source, country"]
```

### 쿼리 작성

```
나: /write-query 가입 월별로 사용자를 그룹화하고 1, 3, 6, 12개월 후에도 여전히 활성
    상태인 비율을 보여주는 코호트 리텐션 분석이 필요합니다. Snowflake를 사용합니다.

Claude: [CTE를 사용해 최적화된 Snowflake SQL 작성]
       → [각 단계에 대한 설명 주석 추가]
       → [파티션 프루닝에 대한 성능 참고 사항 포함]
```

### 대시보드 구축

```
나: /build-dashboard 월간 매출, 주요 제품, 지역별 세부 내역이 포함된 판매 대시보드를
    생성해 주세요. 여기 데이터가 있습니다: [CSV 붙여넣기]

Claude: [독립형 HTML 파일 생성]
       → [대화형 Chart.js 시각화 포함]
       → [지역 및 기간에 대한 드롭다운 필터 추가]
       → [검토를 위해 브라우저에서 열기]
```

### 공유 전 검증

```
나: /validate [분석 문서 공유]

Claude: [방법론 검토] → [이탈 분석에서 생존 편향 여부 확인]
       → [집계 로직 검증] → [플래그: "분모에서 평가판 사용자가 제외되어
          전환율이 약 5pp 높게 측정될 수 있습니다"]
       → [확신도: "유의 사항과 함께 공유 준비 완료"]
```

## 데이터 스택 연결

> 익숙하지 않은 플레이스홀더가 있거나 연결된 도구를 확인해야 하는 경우 [CONNECTORS.md](CONNECTORS.md)를 참조하십시오.

이 플러그인은 데이터 인프라에 연결할 때 가장 효과적입니다. 다음을 위한 MCP 서버를 추가하십시오:

- **데이터 웨어하우스**: Snowflake, Databricks, BigQuery 또는 기타 SQL 호환 데이터베이스
- **분석/BI**: Amplitude, Looker, Tableau 등
- **노트북**: Jupyter, Hex 등
- **스프레드시트**: Google Sheets, Excel
- **데이터 오케스트레이션**: Airflow, dbt, Dagster, Prefect
- **데이터 수집**: Fivetran, Airbyte, Stitch

`.mcp.json` 또는 Claude Code 설정에서 MCP 서버를 구성하여 직접 데이터 접근을 활성화할 수 있습니다.
