---
name: sql-queries
description: 모든 주요 데이터 웨어하우스 방언(Snowflake, BigQuery, Databricks, PostgreSQL 등)에서 정확하고 성능 좋은 SQL을 작성합니다. 쿼리를 작성하거나, 느린 SQL을 최적화하거나, 방언 간 변환을 하거나, CTE, 윈도우 함수, 집계를 활용한 복잡한 분석 쿼리를 구축할 때 사용합니다.
---

# SQL 쿼리 스킬

모든 주요 데이터 웨어하우스 방언에서 정확하고, 성능 좋고, 가독성 높은 SQL을 작성합니다.

## 방언별 참조

### PostgreSQL (Aurora, RDS, Supabase, Neon 포함)

**날짜/시간:**
```sql
-- 현재 날짜/시간
CURRENT_DATE, CURRENT_TIMESTAMP, NOW()

-- 날짜 산술 연산
date_column + INTERVAL '7 days'
date_column - INTERVAL '1 month'

-- 기간 단위로 절삭
DATE_TRUNC('month', created_at)

-- 부분 추출
EXTRACT(YEAR FROM created_at)
EXTRACT(DOW FROM created_at)  -- 0=일요일

-- 형식 지정
TO_CHAR(created_at, 'YYYY-MM-DD')
```

**문자열 함수:**
```sql
-- 결합 (Concatenation)
first_name || ' ' || last_name
CONCAT(first_name, ' ', last_name)

-- 패턴 매칭
column ILIKE '%pattern%'  -- 대소문자 구분 안 함
column ~ '^regex_pattern$'  -- 정규표현식

-- 문자열 조작
LEFT(str, n), RIGHT(str, n)
SPLIT_PART(str, delimiter, position)
REGEXP_REPLACE(str, pattern, replacement)
```

**배열 및 JSON:**
```sql
-- JSON 액세스
data->>'key'  -- 텍스트
data->'nested'->'key'  -- json
data#>>'{path,to,key}'  -- 중첩 텍스트

-- 배열 연산
ARRAY_AGG(column)
ANY(array_column)
array_column @> ARRAY['value']
```

**성능 팁:**
- `EXPLAIN ANALYZE`로 쿼리를 프로파일링합니다
- 자주 필터링/조인하는 컬럼에 인덱스를 생성합니다
- 상관 서브쿼리에는 `IN`보다 `EXISTS`를 사용합니다
- 일반적인 필터 조건에 부분 인덱스를 사용합니다
- 동시 접근에는 커넥션 풀링을 사용합니다

---

### Snowflake

**날짜/시간:**
```sql
-- 현재 날짜/시간
CURRENT_DATE(), CURRENT_TIMESTAMP(), SYSDATE()

-- 날짜 산술 연산
DATEADD(day, 7, date_column)
DATEDIFF(day, start_date, end_date)

-- 기간 단위로 절삭
DATE_TRUNC('month', created_at)

-- 부분 추출
YEAR(created_at), MONTH(created_at), DAY(created_at)
DAYOFWEEK(created_at)

-- 형식 지정
TO_CHAR(created_at, 'YYYY-MM-DD')
```

**문자열 함수:**
```sql
-- 기본적으로 대소문자 구분 안 함 (데이터 정렬에 따라 다름)
column ILIKE '%pattern%'
REGEXP_LIKE(column, 'pattern')

-- JSON 파싱
column:key::string  -- VARIANT를 위한 도트 노테이션
PARSE_JSON('{"key": "value"}')
GET_PATH(variant_col, 'path.to.key')

-- 배열/객체 평면화 (Flatten)
SELECT f.value FROM table, LATERAL FLATTEN(input => array_col) f
```

**반구조화 데이터:**
```sql
-- VARIANT 유형 액세스
data:customer:name::STRING
data:items[0]:price::NUMBER

-- 중첩 구조 평면화
SELECT
    t.id,
    item.value:name::STRING as item_name,
    item.value:qty::NUMBER as quantity
FROM my_table t,
LATERAL FLATTEN(input => t.data:items) item
```

**성능 팁:**
- 대규모 테이블에는 클러스터링 키를 사용합니다 (기존 인덱스가 아닌)
- 파티션 프루닝을 위해 클러스터링 키 컬럼으로 필터링합니다
- 쿼리 복잡도에 맞는 웨어하우스 크기를 설정합니다
- 고비용 쿼리의 재실행을 피하기 위해 `RESULT_SCAN(LAST_QUERY_ID())`을 사용합니다
- 스테이징/임시 데이터에는 transient 테이블을 사용합니다

---

### BigQuery (Google Cloud)

**날짜/시간:**
```sql
-- 현재 날짜/시간
CURRENT_DATE(), CURRENT_TIMESTAMP()

-- 날짜 산술 연산
DATE_ADD(date_column, INTERVAL 7 DAY)
DATE_SUB(date_column, INTERVAL 1 MONTH)
DATE_DIFF(end_date, start_date, DAY)
TIMESTAMP_DIFF(end_ts, start_ts, HOUR)

-- 기간 단위로 절삭
DATE_TRUNC(created_at, MONTH)
TIMESTAMP_TRUNC(created_at, HOUR)

-- 부분 추출
EXTRACT(YEAR FROM created_at)
EXTRACT(DAYOFWEEK FROM created_at)  -- 1=일요일

-- 형식 지정
FORMAT_DATE('%Y-%m-%d', date_column)
FORMAT_TIMESTAMP('%Y-%m-%d %H:%M:%S', ts_column)
```

**문자열 함수:**
```sql
-- ILIKE 없음, LOWER() 사용
LOWER(column) LIKE '%pattern%'
REGEXP_CONTAINS(column, r'pattern')
REGEXP_EXTRACT(column, r'pattern')

-- 문자열 조작
SPLIT(str, delimiter)  -- ARRAY 반환
ARRAY_TO_STRING(array, delimiter)
```

**배열 및 구조체:**
```sql
-- 배열 연산
ARRAY_AGG(column)
UNNEST(array_column)
ARRAY_LENGTH(array_column)
value IN UNNEST(array_column)

-- 구조체(Struct) 액세스
struct_column.field_name
```

**성능 팁:**
- 스캔되는 바이트를 줄이기 위해 항상 파티션 컬럼(보통 날짜)으로 필터링합니다
- 파티션 내에서 자주 필터링하는 컬럼에 클러스터링을 사용합니다
- 대규모 카디널리티 추정에는 `APPROX_COUNT_DISTINCT()`를 사용합니다
- `SELECT *`를 지양합니다 -- 과금이 스캔된 바이트 기준입니다
- 매개변수화된 스크립트에는 `DECLARE`와 `SET`을 사용합니다
- 대규모 쿼리 실행 전에 드라이 런으로 비용을 미리 확인합니다

---

### Redshift (Amazon)

**날짜/시간:**
```sql
-- 현재 날짜/시간
CURRENT_DATE, GETDATE(), SYSDATE

-- 날짜 산술 연산
DATEADD(day, 7, date_column)
DATEDIFF(day, start_date, end_date)

-- 기간 단위로 절삭
DATE_TRUNC('month', created_at)

-- 부분 추출
EXTRACT(YEAR FROM created_at)
DATE_PART('dow', created_at)  -- 0=일요일
```

**문자열 함수:**
```sql
-- 대소문자 구분 안 함
column ILIKE '%pattern%'
REGEXP_INSTR(column, 'pattern') > 0

-- 문자열 조작
SPLIT_PART(str, delimiter, position)
LISTAGG(column, ', ') WITHIN GROUP (ORDER BY column)
```

**성능 팁:**
- 동일 위치 조인을 위해 배포 키를 설계합니다 (DISTKEY)
- 자주 필터링하는 컬럼에 정렬 키를 사용합니다 (SORTKEY)
- `EXPLAIN`으로 쿼리 플랜을 확인합니다
- 크로스 노드 데이터 이동을 지양합니다 (DS_BCAST 및 DS_DIST 주시)
- `ANALYZE`와 `VACUUM`을 정기적으로 실행합니다
- 스키마 유연성을 위해 지연 바인딩 뷰를 사용합니다

---

### Databricks SQL

**날짜/시간:**
```sql
-- 현재 날짜/시간
CURRENT_DATE(), CURRENT_TIMESTAMP()

-- 날짜 산술 연산
DATE_ADD(date_column, 7)
DATEDIFF(end_date, start_date)
ADD_MONTHS(date_column, 1)

-- 기간 단위로 절삭
DATE_TRUNC('MONTH', created_at)
TRUNC(date_column, 'MM')

-- 부분 추출
YEAR(created_at), MONTH(created_at)
DAYOFWEEK(created_at)
```

**Delta Lake 기능:**
```sql
-- Time travel
SELECT * FROM my_table TIMESTAMP AS OF '2024-01-15'
SELECT * FROM my_table VERSION AS OF 42

-- Describe history
DESCRIBE HISTORY my_table

-- Merge (upsert)
MERGE INTO target USING source
ON target.id = source.id
WHEN MATCHED THEN UPDATE SET *
WHEN NOT MATCHED THEN INSERT *
```

**성능 팁:**
- 쿼리 성능을 위해 Delta Lake의 `OPTIMIZE`와 `ZORDER`를 사용합니다
- 컴퓨팅 집약적 쿼리에 Photon 엔진을 활용합니다
- 자주 접근하는 데이터셋에 `CACHE TABLE`을 사용합니다
- 낮은 카디널리티 날짜 컬럼으로 파티셔닝합니다

---

## 일반적인 SQL 패턴

### 윈도우 함수

```sql
-- 순위 지정 (Ranking)
ROW_NUMBER() OVER (PARTITION BY user_id ORDER BY created_at DESC)
RANK() OVER (PARTITION BY category ORDER BY revenue DESC)
DENSE_RANK() OVER (ORDER BY score DESC)

-- 누계 / 이동 평균 (Running totals / moving averages)
SUM(revenue) OVER (ORDER BY date_col ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW) as running_total
AVG(revenue) OVER (ORDER BY date_col ROWS BETWEEN 6 PRECEDING AND CURRENT ROW) as moving_avg_7d

-- 지연(Lag) / 리드(Lead)
LAG(value, 1) OVER (PARTITION BY entity ORDER BY date_col) as prev_value
LEAD(value, 1) OVER (PARTITION BY entity ORDER BY date_col) as next_value

-- 첫 값 / 끝 값 (First / Last value)
FIRST_VALUE(status) OVER (PARTITION BY user_id ORDER BY created_at ROWS BETWEEN UNBOUNDED PRECEDING AND UNBOUNDED FOLLOWING)
LAST_VALUE(status) OVER (PARTITION BY user_id ORDER BY created_at ROWS BETWEEN UNBOUNDED PRECEDING AND UNBOUNDED FOLLOWING)

-- 전체 대비 비율 (Percent of total)
revenue / SUM(revenue) OVER () as pct_of_total
revenue / SUM(revenue) OVER (PARTITION BY category) as pct_of_category
```

### 가독성을 위한 CTE

```sql
WITH
-- 1단계: 기본 모집단 정의
base_users AS (
    SELECT user_id, created_at, plan_type
    FROM users
    WHERE created_at >= DATE '2024-01-01'
      AND status = 'active'
),

-- 2단계: 사용자 수준 지표 계산
user_metrics AS (
    SELECT
        u.user_id,
        u.plan_type,
        COUNT(DISTINCT e.session_id) as session_count,
        SUM(e.revenue) as total_revenue
    FROM base_users u
    LEFT JOIN events e ON u.user_id = e.user_id
    GROUP BY u.user_id, u.plan_type
),

-- 3단계: 요약 수준으로 집계
summary AS (
    SELECT
        plan_type,
        COUNT(*) as user_count,
        AVG(session_count) as avg_sessions,
        SUM(total_revenue) as total_revenue
    FROM user_metrics
    GROUP BY plan_type
)

SELECT * FROM summary ORDER BY total_revenue DESC;
```

### 코호트 리텐션

```sql
WITH cohorts AS (
    SELECT
        user_id,
        DATE_TRUNC('month', first_activity_date) as cohort_month
    FROM users
),
activity AS (
    SELECT
        user_id,
        DATE_TRUNC('month', activity_date) as activity_month
    FROM user_activity
)
SELECT
    c.cohort_month,
    COUNT(DISTINCT c.user_id) as cohort_size,
    COUNT(DISTINCT CASE
        WHEN a.activity_month = c.cohort_month THEN a.user_id
    END) as month_0,
    COUNT(DISTINCT CASE
        WHEN a.activity_month = c.cohort_month + INTERVAL '1 month' THEN a.user_id
    END) as month_1,
    COUNT(DISTINCT CASE
        WHEN a.activity_month = c.cohort_month + INTERVAL '3 months' THEN a.user_id
    END) as month_3
FROM cohorts c
LEFT JOIN activity a ON c.user_id = a.user_id
GROUP BY c.cohort_month
ORDER BY c.cohort_month;
```

### 퍼널 분석

```sql
WITH funnel AS (
    SELECT
        user_id,
        MAX(CASE WHEN event = 'page_view' THEN 1 ELSE 0 END) as step_1_view,
        MAX(CASE WHEN event = 'signup_start' THEN 1 ELSE 0 END) as step_2_start,
        MAX(CASE WHEN event = 'signup_complete' THEN 1 ELSE 0 END) as step_3_complete,
        MAX(CASE WHEN event = 'first_purchase' THEN 1 ELSE 0 END) as step_4_purchase
    FROM events
    WHERE event_date >= CURRENT_DATE - INTERVAL '30 days'
    GROUP BY user_id
)
SELECT
    COUNT(*) as total_users,
    SUM(step_1_view) as viewed,
    SUM(step_2_start) as started_signup,
    SUM(step_3_complete) as completed_signup,
    SUM(step_4_purchase) as purchased,
    ROUND(100.0 * SUM(step_2_start) / NULLIF(SUM(step_1_view), 0), 1) as view_to_start_pct,
    ROUND(100.0 * SUM(step_3_complete) / NULLIF(SUM(step_2_start), 0), 1) as start_to_complete_pct,
    ROUND(100.0 * SUM(step_4_purchase) / NULLIF(SUM(step_3_complete), 0), 1) as complete_to_purchase_pct
FROM funnel;
```

### 중복 제거

```sql
-- Keep the most recent record per key
WITH ranked AS (
    SELECT
        *,
        ROW_NUMBER() OVER (
            PARTITION BY entity_id
            ORDER BY updated_at DESC
        ) as rn
    FROM source_table
)
SELECT * FROM ranked WHERE rn = 1;
```

## 오류 처리 및 디버깅

쿼리가 실패하는 경우:

1. **구문 오류**: 방언별 구문을 확인합니다 (예: BigQuery에서 `ILIKE` 사용 불가, `SAFE_DIVIDE`는 BigQuery에서만)
2. **컬럼을 찾을 수 없음**: 스키마와 대조하여 컬럼명을 확인합니다 -- 오타, 대소문자 구분 확인 (PostgreSQL은 따옴표 붙인 식별자에 대소문자 구분)
3. **유형 불일치**: 다른 유형을 비교할 때 명시적으로 캐스팅합니다 (`CAST(col AS DATE)`, `col::DATE`)
4. **0으로 나누기**: `NULLIF(denominator, 0)` 또는 방언별 안전 나눗셈을 사용합니다
5. **모호한 컬럼**: JOIN에서 항상 테이블 별칭으로 컬럼명을 한정합니다
6. **Group by 오류**: 비집계 컬럼은 모두 GROUP BY에 포함되어야 합니다 (BigQuery는 별칭으로 그룹핑 허용 예외)
