---
description: 생산성 시스템을 초기화하고 대시보드를 엽니다
---

# Start 커맨드

> 낯선 플레이스홀더가 보이거나 어떤 도구가 연결되어 있는지 확인하려면 [CONNECTORS.md](../CONNECTORS.md)를 보세요.

작업과 메모리 시스템을 초기화한 뒤 통합 대시보드를 엽니다.

## 지침

### 1. 존재 여부 확인

작업 디렉터리에서 다음을 확인합니다.
- `TASKS.md` — 작업 목록
- `CLAUDE.md` — 작업 메모리
- `memory/` — 장기 메모리 디렉터리
- `dashboard.html` — 시각 UI

### 2. 누락 항목 생성

**`TASKS.md`가 없으면:** 표준 템플릿( task-management 스킬 참조)으로 생성하고 현재 작업 디렉터리에 둡니다.

**`dashboard.html`이 없으면:** `${CLAUDE_PLUGIN_ROOT}/skills/dashboard.html`에서 현재 작업 디렉터리로 복사합니다.

**`CLAUDE.md`와 `memory/`가 없으면:** 최초 설정입니다. 대시보드를 연 뒤 메모리 부트스트랩 워크플로(아래)를 시작합니다. 현재 작업 디렉터리에 생성합니다.

### 3. 대시보드 열기

`open` 또는 `xdg-open`을 사용하지 마세요. Cowork에서는 에이전트가 VM에서 실행되기 때문에 셸의 open 명령이 사용자 브라우저에 도달하지 않습니다. 대신 다음처럼 안내합니다. "대시보드가 `dashboard.html`에 준비되었습니다. 파일 브라우저에서 열어 시작하세요."

### 4. 사용자 안내

이미 모두 초기화되어 있다면:
```
대시보드를 열었습니다. 작업과 메모리가 모두 로드되었습니다.
- /productivity:update 작업 동기화 및 메모리 확인
- /productivity:update --comprehensive 전체 활동에 대한 심층 스캔
```

메모리가 아직 부트스트랩되지 않았다면 5단계로 진행합니다.

### 5. 메모리 부트스트랩(최초 1회)

`CLAUDE.md`와 `memory/`가 아직 없을 때만 수행합니다.

업무 언어의 최선의 소스는 실제 작업 목록입니다. 실제 업무 = 실제 단축어입니다.

**사용자에게 질문:**
```
할 일 목록은 어디에 있나요? 예:
- 로컬 파일 (예: TASKS.md, todo.txt)
- 앱 (예: Asana, Linear, Jira, Notion, Todoist)
- 노트 파일

작업 목록을 통해 업무 단축어를 학습하겠습니다.
```

**작업 목록에 접근한 후:**

각 작업 항목에서 다음을 점검합니다.
- 별칭일 가능성이 있는 이름
- 약어/두문자
- 프로젝트 참조 또는 코드네임
- 내부 용어/전문 용어

**각 항목은 대화형으로 해석:**

```
Task: "Send PSR to Todd re: Phoenix blockers"

확인하고 싶은 용어가 있습니다:

1. **PSR** - 약어 의미가 무엇인가요?
2. **Todd** - Todd는 누구인가요? (전체 이름, 역할)
3. **Phoenix** - 프로젝트 코드네임인가요? 어떤 프로젝트인가요?
```

이미 해석한 용어는 재질문하지 않고 다음 항목으로 넘어갑니다.

### 6. 선택적 종합 스캔

작업 목록 해석 후 다음을 제안합니다.
```
메시지, 이메일, 문서를 종합 스캔해도 될까요?
시간은 더 걸리지만 사람, 프로젝트, 용어에 대한 더 풍부한 컨텍스트를 쌓을 수 있습니다.

아니면 현재 수준으로 시작하고 나중에 추가해도 됩니다.
```

**종합 스캔을 선택한 경우:**

사용 가능한 MCP 소스에서 데이터를 수집합니다.
- **Chat:** 최근 메시지, 채널, DM
- **Email:** 발신 메시지, 수신자
- **Documents:** 최근 문서, 협업자
- **Calendar:** 미팅, 참석자

발견된 사람/프로젝트/용어를 브레인덤프 형태로 정리하고 신뢰도별로 분류합니다.
- **즉시 추가 가능**(높은 신뢰도) — 바로 추가 제안
- **확인 필요** — 사용자에게 질문
- **빈도 낮음/불명확** — 추후 참고로 보류

### 7. 메모리 파일 작성

수집 결과로 다음을 생성합니다.

**CLAUDE.md** (작업 메모리, 50~80줄):
```markdown
# Memory

## Me
[Name], [Role] on [Team].

## People
| Who | Role |
|-----|------|
| **[Nickname]** | [Full Name], [role] |

## Terms
| Term | Meaning |
|------|---------|
| [acronym] | [expansion] |

## Projects
| Name | What |
|------|------|
| **[Codename]** | [description] |

## Preferences
- [preferences discovered]
```

**memory/** 디렉터리:
- `memory/glossary.md` — 용어/약어/별칭/코드네임 해설
- `memory/people/{name}.md` — 인물 프로필
- `memory/projects/{name}.md` — 프로젝트 상세
- `memory/context/company.md` — 팀, 도구, 프로세스

### 8. 결과 보고

```
생산성 시스템 준비 완료:
- Tasks: TASKS.md (X items)
- Memory: X people, X terms, X projects
- Dashboard: 브라우저에서 열기

/productivity:update로 최신 상태를 유지하세요(--comprehensive로 심층 스캔).
```

## 참고

- 메모리가 이미 초기화되어 있으면 대시보드만 엽니다
- 별칭은 핵심입니다. 사람들이 실제로 어떻게 불리는지 반드시 수집하세요
- 소스를 사용할 수 없으면 건너뛰고 빈 부분을 안내합니다
- 메모리는 부트스트랩 이후에도 자연 대화를 통해 유기적으로 성장합니다
