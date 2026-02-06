---
description: 설정된 템플릿으로 공통 법무 문의 응답 생성
argument-hint: "[inquiry-type]"
---

# /respond -- Generate Response from Templates

> 낯선 플레이스홀더가 보이거나 어떤 도구가 연결되어 있는지 확인하려면 [CONNECTORS.md](../CONNECTORS.md)를 보세요.

템플릿 기반으로 법무 문의 응답 초안을 생성하고, 템플릿 사용이 부적절한 에스컬레이션 조건을 함께 점검합니다.

**중요**: 생성 결과는 법률 자문이 아니며 발송 전 전문가 검토가 필요합니다.

## Invocation

```bash
/respond [inquiry-type]
```

예시 유형:
- `dsr`, `data-subject-request`
- `hold`, `discovery-hold`
- `vendor`, `vendor-question`
- `nda`, `nda-request`
- `privacy`
- `subpoena`
- `insurance`
- `custom`

## Workflow

### 1. 문의 유형 식별

유형이 모호하면 카테고리를 보여주고 확인합니다.

### 2. 템플릿 로드

`legal.local.md` 또는 템플릿 디렉터리에서 템플릿을 찾습니다.
- 있으면 변수와 필수 입력값 식별
- 없으면 기본 구조로 초안 생성 + 템플릿 신설 제안

### 3. 에스컬레이션 트리거 점검

다음이면 템플릿 단독 대응을 피합니다.
- DSR: 미성년자, 규제기관 요청, 보존의무 충돌, 분쟁 당사자
- Hold: 보존 범위 불명확, 형사 리스크 가능성
- Vendor: 소송/해지 위협, 규제준수 이슈, 구속력 있는 약속 위험
- NDA: 경쟁사/정부기밀/M&A 맥락/특수 민감 데이터

탐지 시:
- 왜 위험한지 설명
- 시니어/외부자문 검토 권고
- 최종본 대신 검토용 초안 제공

### 4. 상세 정보 수집

유형별 필수 변수 수집:
- DSR: 요청자, 요청 유형, 적용 규정, 마감
- Hold: 사건명, custodians, 보존 범위, 시행일
- Vendor: 벤더명, 참조 계약, 문의 항목
- NDA: 요청 부서, 상대방, 목적, 상호/편무

### 5. 응답 생성

- 전문적이고 명확한 톤
- 날짜/마감/의무 반영
- 수신자 다음 액션 명확화
- 필요한 주의문구 포함

### 6. 템플릿 생성(없는 경우)

필요 시 템플릿을 새로 작성:
- 카테고리
- Escalation Triggers
- Variables
- Subject
- Body
- Attachments
- Follow-Up

## Output Format

- `To`, `Subject`, 본문
- Escalation Check 결과
- Follow-Up Actions

## Notes

- 발송 전 항상 사용자 리뷰
- 이메일 커넥터가 있으면 초안 메일 생성 제안
- 규제 응답은 마감과 규정 근거를 반드시 표시
