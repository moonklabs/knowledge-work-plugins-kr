---
name: memory-management
description: Claude를 진정한 업무 협력자로 만드는 2단계 기억 시스템입니다. 약어, 두문자어, 별명, 내부 언어를 해독하여 Claude가 동료처럼 요청을 이해합니다. 작업 기억을 위한 CLAUDE.md, 전체 지식 베이스를 위한 memory/ 디렉토리.
---

# 기억 관리

기억은 Claude를 업무 협력자로 만들어줍니다 - 내부 언어를 이해하는 사람입니다.

## 목표

약어를 이해로 변환합니다:

```
사용자: "oracle 건에 대한 PSR을 todd에게 요청해"
              ↓ Claude가 해독
"Todd Martinez(Finance lead)에게 Oracle Systems 거래
 ($2.3M, Q2 마감)에 대한 Pipeline Status Report 준비를 요청합니다"
```

기억이 없으면 이 요청은 의미가 없습니다. 기억이 있으면 Claude는 다음을 알고 있습니다:
- **todd** → Todd Martinez, Finance lead, Slack 선호
- **PSR** → Pipeline Status Report(주간 영업 문서)
- **oracle** → Oracle Systems 거래, 회사가 아님

## 아키텍처

```
CLAUDE.md          ← 핫 캐시(약 30명, 일반 용어)
memory/
  glossary.md      ← 전체 해독 링(모든 것)
  people/          ← 완전한 프로필
  projects/        ← 프로젝트 세부 정보
  context/         ← 회사, 팀, 도구
```

**CLAUDE.md (핫 캐시):**
- 가장 많이 상호작용하는 상위 약 30명
- 약 30개의 가장 일반적인 두문자어/용어
- 활성 프로젝트(5-15개)
- 사용자 선호 사항
- **목표: 일일 해독 필요의 90%를 커버**

**memory/glossary.md (전체 용어집):**
- 완전한 해독 링 - 모든 사람, 모든 용어
- CLAUDE.md에 없는 것을 검색할 때 사용
- 무한정 증가 가능

**memory/people/, projects/, context/:**
- 실행에 필요한 경우 풍부한 세부 정보
- 전체 프로필, 이력, 컨텍스트

## 조회 흐름

```
사용자: "phoenix에 대한 PSR을 todd에게 물어봐"

1. CLAUDE.md 확인(핫 캐시)
   → Todd? ✓ Todd Martinez, Finance
   → PSR? ✓ Pipeline Status Report
   → Phoenix? ✓ DB 마이그레이션 프로젝트

2. 찾을 수 없으면 → memory/glossary.md 검색
   → 전체 용어집에 모든 사람/모든 것 포함

3. 여전히 찾을 수 없으면 → 사용자에게 문의
   → "X가 무엇을 의미하나요? 기억하겠습니다."
```

이 계층화된 접근 방식은 CLAUDE.md를 간결하게 유지(약 100줄)하면서 memory/에서 무제한 확장을 지원합니다.

## 파일 위치

- **작업 기억:** 현재 작업 디렉토리의 `CLAUDE.md`
- **심층 기억:** `memory/` 하위 디렉토리

## 작업 기억 형식 (CLAUDE.md)

압축을 위해 테이블을 사용합니다. 총 약 50-80줄을 목표로 합니다.

```markdown
# Memory

## Me
[이름], [팀]의 [역할]. [내가 하는 일에 대한 한 문장.]

## People
| Who | Role |
|-----|------|
| **Todd** | Todd Martinez, Finance lead |
| **Sarah** | Sarah Chen, Engineering (Platform) |
| **Greg** | Greg Wilson, Sales |
→ 전체 목록: memory/glossary.md, 프로필: memory/people/

## Terms
| Term | Meaning |
|------|---------|
| PSR | Pipeline Status Report |
| P0 | Drop everything priority |
| standup | Daily 9am sync |
→ 전체 용어집: memory/glossary.md

## Projects
| Name | What |
|------|------|
| **Phoenix** | DB migration, Q2 launch |
| **Horizon** | Mobile app redesign |
→ 세부 정보: memory/projects/

## Preferences
- 25분 회의 및 버퍼
- 비동기 우선, 이메일보다 Slack
- 금요일 오후 회의 없음
```

## 심층 기억 형식 (memory/)

**memory/glossary.md** - 해독 링:
```markdown
# Glossary

업무 약어, 두문자어, 내부 언어.

## Acronyms
| Term | Meaning | Context |
|------|---------|---------|
| PSR | Pipeline Status Report | 주간 영업 문서 |
| OKR | Objectives & Key Results | 분기별 계획 |
| P0/P1/P2 | Priority levels | P0 = 모든 것 중단 |

## Internal Terms
| Term | Meaning |
|------|---------|
| standup | #engineering의 일일 9am 동기화 |
| the migration | Project Phoenix 데이터베이스 작업 |
| ship it | 프로덕션에 배포 |
| escalate | 리더십에 루프 |

## Nicknames → Full Names
| Nickname | Person |
|----------|--------|
| Todd | Todd Martinez (Finance) |
| T | Todd Martinez도 |

## Project Codenames
| Codename | Project |
|----------|---------|
| Phoenix | 데이터베이스 마이그레이션 |
| Horizon | 새로운 모바일 앱 |
```

**memory/people/{name}.md:**
```markdown
# Todd Martinez

**Also known as:** Todd, T
**Role:** Finance Lead
**Team:** Finance
**Reports to:** CFO (Michael Chen)

## Communication
- Slack DM 선호
- 빠른 응답, 매우 직접적
- 최적 시간: 아침

## Context
- 모든 PSR 및 재무 보고 처리
- $500k 이상 거래 승인 주요 연락처
- 예측에 대해 Sales와 긴밀히 협력

## Notes
- Cubs 팬, 야구 이야기 좋아함
```

**memory/projects/{name}.md:**
```markdown
# Project Phoenix

**Codename:** Phoenix
**Also called:** "the migration"
**Status:** Active, launching Q2

## What It Is
레거시 Oracle에서 PostgreSQL로의 데이터베이스 마이그레이션.

## Key People
- Sarah - 기술 리드
- Todd - 예산 소유자
- Greg - 이해관계자(영업 영향)

## Context
$1.2M 예산, 6개월 일정. Horizon 프로젝트의 중요 경로.
```

**memory/context/company.md:**
```markdown
# Company Context

## Tools & Systems
| Tool | Used for | Internal name |
|------|----------|---------------|
| Slack | Communication | - |
| Asana | Engineering tasks | - |
| Salesforce | CRM | "SF" 또는 "the CRM" |
| Notion | Docs/wiki | - |

## Teams
| Team | What they do | Key people |
|------|--------------|------------|
| Platform | Infrastructure | Sarah (lead) |
| Finance | Money stuff | Todd (lead) |
| Sales | Revenue | Greg |

## Processes
| Process | What it means |
|---------|---------------|
| Weekly sync | 월요일 10am 전체 회의 |
| Ship review | 목요일 배포 승인 |
```

## 상호작용 방법

### 사용자 입력 해독 (계층화된 조회)

요청에 따라 행동하기 전에 **항상** 약어를 해독합니다:

```
1. CLAUDE.md (핫 캐시)     → 먼저 확인, 케이스의 90% 커버
2. memory/glossary.md        → 핫 캐시에 없으면 전체 용어집
3. memory/people/, projects/ → 필요할 때 풍부한 세부 정보
4. 사용자에게 문의          → 알 수 없는 용어? 배웁니다.
```

예시:
```
사용자: "oracle에 대한 PSR을 todd에게 요청해"

CLAUDE.md 조회:
  "todd" → Todd Martinez, Finance ✓
  "PSR" → Pipeline Status Report ✓
  "oracle" → (핫 캐시에 없음)

memory/glossary.md 조회:
  "oracle" → Oracle Systems 거래 ($2.3M) ✓

이제 Claude는 전체 컨텍스트로 행동할 수 있습니다.
```

### 기억 추가

사용자가 "이걸 기억해" 또는 "X는 Y를 의미해"라고 말할 때:

1. **용어집 항목**(두문자어, 용어, 약어):
   - memory/glossary.md에 추가
   - 자주 사용되면 CLAUDE.md Quick Glossary에 추가

2. **사람:**
   - memory/people/{name}.md 생성/업데이트
   - 중요하면 CLAUDE.md Key People에 추가
   - **별명 캡처** - 해독에 중요

3. **프로젝트:**
   - memory/projects/{name}.md 생성/업데이트
   - 현재 진행 중이면 CLAUDE.md Active Projects에 추가
   - **코드명 캡처** - "Phoenix", "the migration" 등

4. **선호 사항:** CLAUDE.md Preferences 섹션에 추가

### 기억 회상

사용자가 "X가 누구야" 또는 "X가 무엇을 의미해"라고 물을 때:

1. 먼저 CLAUDE.md 확인
2. 전체 세부 정보는 memory/ 확인
3. 찾을 수 없으면: "X가 무엇을 의미하는지 아직 모릅니다. 알려주시겠어요?"

### 점진적 공개

1. 요청의 빠른 파싱을 위해 CLAUDE.md 로드
2. 실행을 위한 전체 컨텍스트가 필요할 때 memory/에 진입
3. 예시: todd에게 PSR에 대한 이메일 초안 작성
   - CLAUDE.md는 Todd = Todd Martinez, PSR = Pipeline Status Report를 알려줌
   - memory/people/todd-martinez.md는 그가 Slack을 선호하고 직접적임을 알려줌

## 부트스트랩

`/productivity:start`를 사용하여 채팅, 캘린더, 이메일, 문서를 스캔하여 초기화합니다. 사람, 프로젝트를 추출하고 용어집 작성을 시작합니다.

## 규칙

- 가독성을 위해 CLAUDE.md의 용어를 **굵게** 표시
- CLAUDE.md를 약 100줄 이하로 유지("핫 30" 규칙)
- 파일명: 소문자, 하이픈(`todd-martinez.md`, `project-phoenix.md`)
- 항상 별명과 대체 이름 캡처
- 쉬운 조회를 위한 용어집 테이블
- 자주 사용되면 CLAUDE.md로 승격
- 오래되면 memory/로만 강등

## 어디에 무엇을 넣을지

| 유형 | CLAUDE.md (핫 캐시) | memory/ (전체 저장소) |
|------|---------------------|---------------------|
| 사람 | 상위 약 30명의 빈번한 연락처 | glossary.md + people/{name}.md |
| 두문자어/용어 | 약 30개의 가장 일반적인 것 | glossary.md (완전한 목록) |
| 프로젝트 | 활성 프로젝트만 | glossary.md + projects/{name}.md |
| 별명 | 상위 30에 있으면 Key People에 | glossary.md (모든 별명) |
| 회사 컨텍스트 | 빠른 참조만 | context/company.md |
| 선호 사항 | 모든 선호 사항 | - |
| 이력/오래된 것 | ✗ 제거 | ✓ memory/에 보관 |

## 승격 / 강등

**CLAUDE.md로 승격할 때:**
- 용어/사람을 자주 사용
- 활성 작업의 일부

**memory/로만 강등할 때:**
- 프로젝트 완료
- 더 이상 빈번한 연락처가 아닌 사람
- 거의 사용되지 않는 용어

이렇게 하면 CLAUDE.md가 신선하고 관련성 있게 유지됩니다.
