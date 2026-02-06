---
description: 작업을 동기화하고 현재 활동에서 메모리를 갱신
argument-hint: "[--comprehensive]"
---

# Update 커맨드

> 낯선 플레이스홀더가 보이거나 어떤 도구가 연결되어 있는지 확인하려면 [CONNECTORS.md](../CONNECTORS.md)를 보세요.

작업 목록과 메모리를 최신 상태로 유지합니다. 두 가지 모드가 있습니다.

- **기본:** 외부 도구에서 작업을 동기화하고, 오래된 항목을 분류하고, 메모리 빈틈을 점검
- **`--comprehensive`:** 채팅/이메일/캘린더/문서를 심층 스캔해 누락된 TODO를 표시하고 새 메모리를 제안

## 사용법

```bash
/productivity:update
/productivity:update --comprehensive
```

## 기본 모드

### 1. 현재 상태 로드

`TASKS.md`와 `memory/` 디렉터리를 읽습니다. 없으면 `/productivity:start`를 먼저 권장합니다.

### 2. 외부 소스에서 작업 동기화

사용 가능한 작업 소스를 확인합니다.
- **프로젝트 트래커**(예: Asana, Linear, Jira) (MCP 사용 가능 시)
- **GitHub Issues** (레포 내): `gh issue list --assignee=@me`

사용 가능한 소스가 없으면 3단계로 넘어갑니다.

**사용자에게 할당된 작업**(open/in-progress)을 가져와 TASKS.md와 비교합니다.

| 외부 작업 | TASKS.md 매칭? | 조치 |
|---------------|-----------------|--------|
| 외부에 있고 TASKS.md 없음 | 매칭 없음 | 추가 제안 |
| 외부에 있고 TASKS.md에 있음 | 제목 유사 매칭 | 건너뜀 |
| TASKS.md에 있고 외부에 없음 | 매칭 없음 | 오래된 항목으로 플래그 |
| 외부에서 완료됨 | Active 섹션에 있음 | 완료 처리 제안 |

차이를 보여주고 무엇을 추가/완료할지 사용자가 결정하도록 합니다.

### 3. 오래된 항목 분류

TASKS.md의 Active 작업을 검토하고 다음을 표시합니다.
- 마감일이 지난 작업
- Active 상태로 30일 이상인 작업
- 컨텍스트 없는 작업(사람/프로젝트 없음)

각 항목에 대해 분류를 요청합니다: 완료 처리? 재일정? Someday로 이동?

### 4. 메모리 빈틈을 위한 작업 해석

각 작업에서 모든 엔터티(사람, 프로젝트, 약어, 도구, 링크)를 해석합니다.

```
Task: "Send PSR to Todd re: Phoenix blockers"

Decode:
- PSR → ✓ Pipeline Status Report (in glossary)
- Todd → ✓ Todd Martinez (in people/)
- Phoenix → ? Not in memory
```

완전히 해석된 항목과 빈틈이 있는 항목을 추적합니다.

### 5. 빈틈 채우기

알 수 없는 용어를 묶어서 제시합니다.
```
작업에 있는 용어 중 컨텍스트가 없는 항목이 있습니다:

1. "Phoenix" (from: "Send PSR to Todd re: Phoenix blockers")
   → Phoenix가 무엇인가요?

2. "Maya" (from: "sync with Maya on API design")
   → Maya는 누구인가요?
```

답변을 적절한 메모리 파일(people/, projects/, glossary.md)에 추가합니다.

### 6. 컨텍스트 보강

작업에는 메모리보다 풍부한 컨텍스트가 포함됩니다. 이를 추출해 업데이트합니다.
- **링크** → 프로젝트/인물 파일에 추가
- **상태 변경**("launch done") → 프로젝트 상태 업데이트, CLAUDE.md에서 비중 낮추기
- **관계**("Todd's sign-off on Maya's proposal") → 인물 간 참조 추가
- **마감** → 프로젝트 파일에 추가

### 7. 보고

```
업데이트 완료:
- Tasks: 프로젝트 트래커(예: Asana)에서 +3, 1개 완료, 2개 분류
- Memory: 2개 빈틈 보완, 1개 프로젝트 보강
- 모든 작업 해석 완료 ✓
```

## 종합 모드 (`--comprehensive`)

기본 모드의 모든 작업 + 최근 활동에 대한 심층 스캔을 수행합니다.

### 추가 단계: 활동 소스 스캔

사용 가능한 MCP 소스에서 데이터를 수집합니다.
- **Chat:** 최근 메시지 검색, 활성 채널 읽기
- **Email:** 보낸 메일 검색
- **Documents:** 최근 변경 문서 목록
- **Calendar:** 최근 + 예정된 이벤트 목록

### 추가 단계: 누락된 TODO 표시

활동을 TASKS.md와 비교해 누락된 액션 아이템을 노출합니다.

```
## 누락 가능 작업

활동에서 아직 캡처되지 않은 TODO로 보이는 항목:

1. 채팅에서 (Jan 18):
   "I'll send the updated mockups by Friday"
   → TASKS.md에 추가할까요?

2. "Phoenix Standup" 미팅 (Jan 17):
   반복 미팅이 있지만 Phoenix 작업이 없음
   → 필요한 작업이 있나요?

3. 이메일에서 (Jan 16):
   "I'll review the API spec this week"
   → TASKS.md에 추가할까요?
```

사용자가 추가할 항목을 선택합니다.

### 추가 단계: 새 메모리 제안

메모리에 없는 새 엔터티를 노출합니다.

```
## 새 인물(메모리 없음)
| Name | Frequency | Context |
|------|-----------|---------|
| Maya Rodriguez | 12 mentions | design, UI reviews |
| Alex K | 8 mentions | DMs about API |

## 새 프로젝트/주제
| Name | Frequency | Context |
|------|-----------|---------|
| Starlight | 15 mentions | planning docs, product |

## 정리 제안
- **Horizon project** — 30일 동안 언급 없음. 완료로 표시할까요?
```

신뢰도별로 묶어 제시합니다. 높은 신뢰도는 바로 추가를 제안하고, 낮은 신뢰도는 확인 질문을 합니다.

## 참고

- 사용자 확인 없이 작업이나 메모리를 자동으로 추가하지 않습니다
- 외부 소스 링크는 가능한 경우 유지합니다
- 작업 제목은 유사도 매칭으로 경미한 문구 차이를 처리합니다
- 자주 실행해도 안전합니다 — 새로운 정보가 있을 때만 업데이트됩니다
- `--comprehensive`는 항상 대화형으로 진행됩니다
