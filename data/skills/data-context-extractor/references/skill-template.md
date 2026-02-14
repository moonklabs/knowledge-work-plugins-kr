# 생성된 스킬 템플릿

새로운 데이터 분석 스킬을 생성할 때 이 템플릿을 사용합니다. 모든 `[PLACEHOLDER]` 값을 교체하십시오.

---

```markdown
---
name: [company]-data-analyst
description: "[COMPANY] 데이터 분석 스킬. 엔티티 정의, 지표 계산 및 일반적인 쿼리 패턴을 포함하여 [WAREHOUSE_TYPE] 쿼리를 위한 컨텍스트를 제공합니다. 다음 항목을 위해 [COMPANY] 데이터를 분석할 때 사용하십시오: (1) [PRIMARY_USE_CASE_1], (2) [PRIMARY_USE_CASE_2], (3) [PRIMARY_USE_CASE_3] 또는 [COMPANY] 고유의 컨텍스트가 필요한 모든 데이터 질문."
---

# [COMPANY] 데이터 분석

## SQL 방언: [WAREHOUSE_TYPE]

[sql-dialects.md에서 적절한 방언 섹션을 삽입]

---

## 엔티티 구분

사용자가 다음 용어를 언급할 때 어떤 엔티티를 의미하는지 명확히 합니다:

[예시 형식 - 탐색 결과에 따라 커스터마이징:]

**"사용자"의 의미:**
- **계정(Account)**: 개별 로그인/프로필 ([PRIMARY_TABLE]: [ID_FIELD])
- **조직(Organization)**: 여러 계정을 가질 수 있는 청구 엔티티 ([ORG_TABLE]: [ORG_ID])
- **[OTHER_TYPE]**: [DEFINITION] ([TABLE]: [ID])

**관계:**
- [ENTITY_1] → [ENTITY_2]: [RELATIONSHIP_TYPE] ([JOIN_KEY]로 조인)

---

## 비즈니스 용어

| 용어 | 정의 | 참고 |
|------|------|------|
| [TERM_1] | [DEFINITION] | [CONTEXT/GOTCHA] |
| [TERM_2] | [DEFINITION] | [CONTEXT/GOTCHA] |
| [ACRONYM] | [FULL_NAME] - [EXPLANATION] | |

---

## 표준 필터

명시적으로 달리 지시하지 않는 한 항상 다음 필터를 적용합니다:

```sql
-- 테스트/내부 데이터 제외
WHERE [TEST_FLAG_COLUMN] = FALSE
  AND [INTERNAL_FLAG_COLUMN] = FALSE

-- 유효하지 않은/사기 데이터 제외
  AND [STATUS_COLUMN] != '[EXCLUDED_STATUS]'

-- [기타 표준 제외 사항]
```

**재정의 시점:**
- [SCENARIO_1]: [CONDITION]일 때 [NORMALLY_EXCLUDED]를 포함

---

## 핵심 지표

### [METRIC_1_NAME]
- **정의**: [PLAIN_ENGLISH_EXPLANATION]
- **공식**: `[EXACT_CALCULATION]`
- **소스**: `[TABLE_NAME].[COLUMN_NAME]`
- **시간 단위**: [DAILY/WEEKLY/MONTHLY]
- **주의사항**: [EDGE_CASES_OR_GOTCHAS]

### [METRIC_2_NAME]
[형식 반복]

---

## 데이터 최신성

| 테이블 | 업데이트 빈도 | 일반적인 지연 |
|--------|-------------|-------------|
| [TABLE_1] | [FREQUENCY] | [LAG] |
| [TABLE_2] | [FREQUENCY] | [LAG] |

데이터 최신성 확인:
```sql
SELECT MAX([DATE_COLUMN]) as latest_data FROM [TABLE]
```

---

## 지식 베이스 탐색

상세 테이블 문서화를 위해 다음 참조 파일을 사용합니다:

| 도메인 | 참조 파일 | 용도 |
|--------|----------|------|
| [DOMAIN_1] | `references/[domain1].md` | [BRIEF_DESCRIPTION] |
| [DOMAIN_2] | `references/[domain2].md` | [BRIEF_DESCRIPTION] |
| 엔티티 | `references/entities.md` | 엔티티 정의 및 관계 |
| 지표 | `references/metrics.md` | KPI 계산 및 공식 |

---

## 일반적인 쿼리 패턴

### [PATTERN_1_NAME]
```sql
[SAMPLE_QUERY]
```

### [PATTERN_2_NAME]
```sql
[SAMPLE_QUERY]
```

---

## 문제 해결

### 일반적인 실수
- **[MISTAKE_1]**: [EXPLANATION] → [CORRECT_APPROACH]
- **[MISTAKE_2]**: [EXPLANATION] → [CORRECT_APPROACH]

### 접근 권한 문제
- `[TABLE]`에서 권한 오류가 발생하면: [WORKAROUND]
- PII 제한 컬럼의 경우: [ALTERNATIVE_APPROACH]

### 성능 팁
- 스캔되는 데이터를 줄이기 위해 먼저 `[PARTITION_COLUMN]`으로 필터링합니다
- 대규모 테이블의 경우 탐색 시 `LIMIT`을 사용합니다
- 가능한 경우 `[RAW_TABLE]`보다 `[AGGREGATED_TABLE]`을 선호합니다
```

---

## 커스터마이징 참고사항

스킬 생성 시:

1. **모든 플레이스홀더를 채우십시오** - `[PLACEHOLDER]` 텍스트를 남기지 마십시오
2. **사용하지 않는 섹션을 제거하십시오** - 대시보드가 없으면 해당 섹션을 제거합니다
3. **구체성을 추가하십시오** - 일반적인 조언보다 구체적인 컬럼명과 값이 더 유용합니다
4. **실제 예시를 포함하십시오** - 샘플 쿼리는 실제 테이블/컬럼명을 사용해야 합니다
5. **스캔 가능하게 유지하십시오** - 테이블과 코드 블록을 적극적으로 활용합니다
