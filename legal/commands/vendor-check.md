---
description: 연결된 시스템 전반에서 특정 벤더와의 기존 계약 상태 점검
argument-hint: "[vendor name]"
---

# /vendor-check -- Vendor Agreement Status

> 낯선 플레이스홀더가 보이거나 어떤 도구가 연결되어 있는지 확인하려면 [CONNECTORS.md](../CONNECTORS.md)를 보세요.

연결된 시스템에서 벤더 관련 계약 현황을 통합 조회합니다.

**중요**: 본 결과는 참고용이며 원문 계약서를 기준으로 전문가 검증이 필요합니다.

## Invocation

```bash
/vendor-check [vendor name]
```

## Workflow

### 1. 벤더 식별

- 법인명/상호/약칭/모회사-자회사 관계를 정규화
- 모호하면 사용자 확인

### 2. 시스템별 검색

우선순위:
1. CLM: 체결/만료/협상중/개정
2. CRM: 계정 상태, 딜, 담당자
3. Email: 최근 협상/첨부/스레드
4. Documents: 체결본/레드라인/실사자료
5. Chat: 최근 요청/질문/논의

### 3. 계약 상태 컴파일

계약별로:
- Agreement Type
- Status
- Effective Date
- Expiration Date
- Auto-Renewal
- Key Terms
- Amendments

### 4. 갭 분석

체크리스트:
- NDA
- MSA
- DPA
- SOW
- SLA
- 보험증권

관계 유형 대비 필요한 문서가 누락되면 플래그합니다.

### 5. 리포트 출력

- 검색일
- 확인한 소스/미연결 소스
- 관계 개요
- 계약 요약
- 갭 분석
- 임박 만료/갱신 액션
- 이메일/채팅 컨텍스트 메모

### 6. 미연결 소스 처리

핵심 소스가 없으면 누락 범위를 명확히 표시하고 수동 확인 경로를 안내합니다.

## Notes

- 어떤 소스에서도 계약이 없으면 명확히 보고하고 다른 저장 위치 확인 요청
- 90일 이내 만료는 상단 강조
- 만료 계약의 생존 의무(비밀유지/면책 등) 여부를 별도 표기
