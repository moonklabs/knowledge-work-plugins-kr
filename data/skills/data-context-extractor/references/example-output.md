# 예시: 생성된 스킬

부트스트랩 프로세스 후 생성된 스킬의 예시입니다. 이 예시는 Snowflake를 사용하는 가상의 전자상거래 회사 "ShopCo"를 위한 것입니다.

---

## 예시 SKILL.md

```markdown
---
name: shopco-data-analyst
description: "Snowflake용 ShopCo 데이터 분석 스킬. 고객, 주문 및 제품 분석을 포함한 전자상거래 데이터 쿼리에 대한 컨텍스트를 제공합니다. 다음 항목을 위해 ShopCo 데이터를 분석할 때 사용하십시오: (1) 수익 및 주문 지표, (2) 고객 행동 및 리텐션, (3) 제품 성능 또는 ShopCo 고유의 컨텍스트가 필요한 모든 데이터 질문."
---

# ShopCo 데이터 분석

## SQL 방언: Snowflake

- **테이블 참조**: `SHOPCO_DW.SCHEMA.TABLE` 또는 대소문자 구분 시 따옴표 사용: `"Column_Name"`
- **안전한 나눗셈**: `DIV0(a, b)`는 0 반환, `DIV0NULL(a, b)`는 NULL 반환
- **날짜 함수**:
  - `DATE_TRUNC('MONTH', date_col)`
  - `DATEADD(DAY, -1, date_col)`
  - `DATEDIFF(DAY, start_date, end_date)`
- **컬럼 제외**: `SELECT * EXCLUDE (column_to_exclude)`

---

## 엔티티 구분

**"고객"의 의미:**
- **사용자(User)**: 탐색 및 아이템 저장이 가능한 로그인 계정 (CORE.DIM_USERS: user_id)
- **고객(Customer)**: 최소 한 번 이상 구매한 사용자 (CORE.DIM_CUSTOMERS: customer_id)
- **계정(Account)**: 청구 엔티티로, B2B에서 여러 사용자를 가질 수 있음 (CORE.DIM_ACCOUNTS: account_id)

**관계:**
- 사용자 → 고객: 1:1 (구매자의 경우 customer_id = user_id)
- 계정 → 사용자: 1:다 (account_id로 조인)

---

## 비즈니스 용어

| 용어 | 정의 | 참고 |
|------|------|------|
| GMV | 총 상품 거래액(Gross Merchandise Value) - 반품/할인 전 총 주문 금액 | 톱라인 보고에 사용 |
| NMV | 순 상품 거래액(Net Merchandise Value) - GMV에서 반품과 할인을 차감 | 실제 매출에 사용 |
| AOV | 평균 주문 금액(Average Order Value) - NMV / 주문 수 | $0 주문 제외 |
| LTV | 생애 가치(Lifetime Value) - 첫 주문 이후 고객당 총 NMV | 롤링 계산, 일별 업데이트 |
| CAC | 고객 획득 비용(Customer Acquisition Cost) - 마케팅 비용 / 신규 고객 | 코호트 월별 |

---

## 표준 필터

명시적으로 달리 지시하지 않는 한 항상 다음 필터를 적용합니다:

```sql
-- 테스트 및 내부 주문 제외
WHERE order_status != 'TEST'
  AND customer_type != 'INTERNAL'
  AND is_employee_order = FALSE

-- 수익 지표에서 취소된 주문 제외
  AND order_status NOT IN ('CANCELLED', 'FRAUDULENT')
```

---

## 핵심 지표

### 총 상품 거래액 (GMV)
- **정의**: 모든 주문의 총 금액
- **공식**: `SUM(order_total_gross)`
- **소스**: `CORE.FCT_ORDERS.order_total_gross`
- **시간 단위**: 일별, 주별/월별로 집계
- **주의사항**: 이후 취소 또는 반품될 수 있는 주문을 포함

### 순 매출
- **정의**: 반품 및 할인 후 실제 매출
- **공식**: `SUM(order_total_gross - return_amount - discount_amount)`
- **소스**: `CORE.FCT_ORDERS`
- **주의사항**: 반품은 주문 후 최대 90일까지 발생 가능; 확정 수치는 settled_revenue 사용

---

## 지식 베이스 탐색

| 도메인 | 참조 파일 | 용도 |
|--------|----------|------|
| 주문 | `references/orders.md` | 주문 테이블, GMV/NMV 계산 |
| 고객 | `references/customers.md` | 사용자/고객 엔티티, LTV, 코호트 |
| 제품 | `references/products.md` | 카탈로그, 재고, 카테고리 |

---

## 일반적인 쿼리 패턴

### 채널별 일별 GMV
```sql
SELECT
    DATE_TRUNC('DAY', order_timestamp) AS order_date,
    channel,
    SUM(order_total_gross) AS gmv,
    COUNT(DISTINCT order_id) AS order_count
FROM SHOPCO_DW.CORE.FCT_ORDERS
WHERE order_status NOT IN ('TEST', 'CANCELLED', 'FRAUDULENT')
  AND order_timestamp >= DATEADD(DAY, -30, CURRENT_DATE())
GROUP BY 1, 2
ORDER BY 1 DESC, 3 DESC
```

### 고객 코호트 리텐션
```sql
WITH cohorts AS (
    SELECT
        customer_id,
        DATE_TRUNC('MONTH', first_order_date) AS cohort_month
    FROM SHOPCO_DW.CORE.DIM_CUSTOMERS
)
SELECT
    c.cohort_month,
    DATEDIFF(MONTH, c.cohort_month, DATE_TRUNC('MONTH', o.order_timestamp)) AS months_since_first,
    COUNT(DISTINCT c.customer_id) AS active_customers
FROM cohorts c
JOIN SHOPCO_DW.CORE.FCT_ORDERS o ON c.customer_id = o.customer_id
WHERE o.order_status NOT IN ('TEST', 'CANCELLED')
GROUP BY 1, 2
ORDER BY 1, 2
```
```

---

## 예시 references/orders.md

```markdown
# 주문 테이블

ShopCo의 주문 및 거래 데이터입니다.

---

## 핵심 테이블

### FCT_ORDERS
**위치**: `SHOPCO_DW.CORE.FCT_ORDERS`
**설명**: 모든 주문의 팩트 테이블. 주문당 하나의 행.
**기본 키**: `order_id`
**업데이트 빈도**: 시간별 (15분 지연)
**파티션 기준**: `order_date`

| 컬럼 | 유형 | 설명 | 참고 |
|------|------|------|------|
| **order_id** | VARCHAR | 고유 주문 식별자 | |
| **customer_id** | VARCHAR | DIM_CUSTOMERS FK | 비회원 결제 시 NULL |
| **order_timestamp** | TIMESTAMP_NTZ | 주문 시점 | UTC |
| **order_date** | DATE | order_timestamp의 날짜 부분 | 파티션 컬럼 |
| **order_status** | VARCHAR | 현재 상태 | PENDING, SHIPPED, DELIVERED, CANCELLED, RETURNED |
| **channel** | VARCHAR | 유입 채널 | WEB, APP, MARKETPLACE |
| **order_total_gross** | DECIMAL(12,2) | 할인 전 총액 | |
| **discount_amount** | DECIMAL(12,2) | 적용된 총 할인 | |
| **return_amount** | DECIMAL(12,2) | 반품된 품목 금액 | 비동기 업데이트 |

**관계**:
- `customer_id`로 `DIM_CUSTOMERS`와 조인
- `order_id`를 통해 `FCT_ORDER_ITEMS`의 상위

---

## 샘플 쿼리

### 반품률이 포함된 주문
```sql
SELECT
    DATE_TRUNC('WEEK', order_date) AS week,
    COUNT(*) AS total_orders,
    SUM(CASE WHEN return_amount > 0 THEN 1 ELSE 0 END) AS orders_with_returns,
    DIV0(SUM(CASE WHEN return_amount > 0 THEN 1 ELSE 0 END), COUNT(*)) AS return_rate
FROM SHOPCO_DW.CORE.FCT_ORDERS
WHERE order_status NOT IN ('TEST', 'CANCELLED')
  AND order_date >= DATEADD(MONTH, -3, CURRENT_DATE())
GROUP BY 1
ORDER BY 1
```
```

---

이 예시는 다음을 보여줍니다:
- 트리거 설명이 포함된 완전한 프론트매터
- 방언별 SQL 참고사항
- 명확한 엔티티 구분
- 용어 사전
- 복사하여 사용할 수 있는 SQL 형태의 표준 필터
- 공식이 포함된 지표 정의
- 참조 파일로의 탐색
- 실제 실행 가능한 쿼리 예시
