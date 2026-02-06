---
description: 법무 업무용 컨텍스트 브리핑 생성(데일리/주제/인시던트)
argument-hint: "[daily | topic <query> | incident]"
---

# /brief -- Legal Team Briefing

> 낯선 플레이스홀더가 보이거나 어떤 도구가 연결되어 있는지 확인하려면 [CONNECTORS.md](../CONNECTORS.md)를 보세요.

법무 업무 브리핑을 생성합니다. `daily`, `topic`, `incident` 모드를 지원합니다.

**중요**: 이 커맨드는 법무 워크플로를 보조하지만 법률 자문을 제공하지 않습니다. 결과는 법무 전문가 검토 후 활용해야 합니다.

## Invocation

```bash
/brief daily
/brief topic [query]
/brief incident [topic]
```

모드가 없으면 필요한 브리프 유형을 확인합니다.

## Modes

### Daily Brief

아침 업무 시작용 요약:
- 이메일: 신규 계약요청, 컴플라이언스 문의, 협상 회신, 외부자문 메일
- 캘린더: 법무 준비가 필요한 회의, 임박 마감
- 채팅: 법무 관련 멘션/DM/에스컬레이션
- CLM: 검토 대기, 서명 대기, 만료 임박 계약
- CRM: 법무 개입이 필요한 단계로 이동한 딜

출력:
- 긴급/즉시 조치
- 계약 파이프라인
- 신규 요청
- 오늘 일정
- 팀 활동
- 주간 마감
- 미연결/오류 소스

### Topic Brief

특정 법무 주제 조사 브리프:
1. 사용자 쿼리 수집
2. 문서/이메일/채팅/CLM 검색
3. 요약 생성

출력:
- Summary
- Background
- Current State
- Key Considerations
- Internal Precedent
- Gaps
- Recommended Next Steps

주의:
- 내부 연결 소스 기반 요약이며 정식 법률조사 대체 아님
- 판례/최신 법령이 필요하면 전문 리서치 플랫폼 또는 외부자문 권고

### Incident Brief

데이터 유출, 소송 위협, 규제 질의 등 긴급 사안 대응 브리프:
1. 사안 설명 수집
2. 이메일/채팅/문서/캘린더/CLM 신속 스캔
3. 즉시 실행 가능한 브리프 생성

출력:
- Situation Summary
- Timeline
- Immediate Legal Considerations
- Relevant Agreements
- Internal Response
- Key Contacts
- Recommended Immediate Actions
- Information Gaps
- Sources Checked

인시던트 브리프 원칙:
- 완전성보다 속도 우선
- 보존의무/리티게이션 홀드 즉시 플래그
- privilege 고려사항 명시
- 데이터 유출 가능 시 통지 기한 플래그

## General Notes

- 미연결 소스는 눈에 띄게 표기
- 항상 실행 가능한 다음 단계 중심으로 작성
- 길게 인용하지 말고 핵심과 링크 중심으로 요약
