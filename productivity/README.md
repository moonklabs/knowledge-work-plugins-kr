# Productivity 플러그인

생산성 플러그인입니다. 주로 Anthropic의 agentic 데스크톱 애플리케이션인 [Cowork](https://claude.com/product/cowork)용으로 설계되었지만 Claude Code에서도 동작합니다. 작업 관리, 업무 메모리, 시각적 대시보드를 제공하며, Claude가 사람, 프로젝트, 용어를 학습해 chatbot이 아니라 동료처럼 일할 수 있게 합니다.

## 설치

```
claude plugins add knowledge-work-plugins/productivity
```

## 주요 기능

이 플러그인은 Claude가 사용자의 일을 지속적으로 이해하도록 돕습니다.

- **작업 관리** — Claude가 읽고 쓰며 실행 기준으로 삼는 markdown 작업 목록(`TASKS.md`)입니다. 자연스럽게 작업을 추가하면 Claude가 상태를 추적하고 오래된 항목을 분류하며 외부 도구와 동기화합니다.
- **업무 메모리** — 약어, 사람, 프로젝트, 용어를 Claude에게 가르치는 2계층 메모리 시스템입니다. "Todd에게 Oracle PSR 해달라고 해"라고 말하면 Claude가 누구에게 무엇을 어떤 거래에 대해 요청하는지 이해합니다.
- **시각적 대시보드** — 작업 보드 보기와 Claude가 workplace에 대해 아는 내용을 보여주는 local HTML file입니다. 보드나 파일 어느 쪽에서 수정해도 동기화됩니다.

## 명령

| 명령 | 설명 |
|---------|--------------|
| `/start` | 작업과 메모리를 초기화하고 dashboard를 엽니다. |
| `/update` | 오래된 항목을 분류하고 memory gap을 확인하며, 가능한 경우 외부 도구에서 동기화합니다. |
| `/update --comprehensive` | Email, calendar, chat을 심층 스캔해 놓친 할 일을 표시하고 새 메모리를 제안합니다. |

## 스킬

| 스킬 | 설명 |
|-------|-------------|
| `memory-management` | Claude를 실제 업무 협업자처럼 만드는 2계층 메모리 시스템입니다. 약어, 두문자어, 별칭, 내부 언어를 해석해 Claude가 동료처럼 요청을 이해하게 합니다. 작업 메모리에는 CLAUDE.md를, 전체 지식 기반에는 memory/ 디렉터리를 사용합니다. |
| `task-management` | 공유 TASKS.md 파일을 사용하는 간단한 작업 관리입니다. 사용자가 자신의 작업을 묻거나, 작업 추가/완료를 원하거나, 약속 추적에 도움이 필요할 때 참조합니다. |

## 예시 워크플로

### 시작하기

```
You: /start

Claude: [TASKS.md, CLAUDE.md, memory/ directory, dashboard.html을 생성]
        [browser에서 dashboard를 엽니다]
        [memory seed를 위해 역할, 팀, 현재 우선순위를 질문]
```

### 자연어로 작업 추가

```
You: 금요일까지 Sarah의 budget proposal을 review해야 하고,
     Greg와 sync한 뒤 Q2 roadmap을 draft해야 하고,
     Platform team의 API spec도 follow up해야 해

Claude: [세 task를 context와 함께 TASKS.md에 추가]
        [Dashboard가 자동으로 업데이트]
```

### 아침 동기화

```
You: /update --comprehensive

Claude: [email, calendar, chat에서 새 action item을 scan]
        [Flags: "Budget proposal review is due tomorrow — still open"]
        [Suggests: "New person mentioned in 3 threads: Jamie Park,
         Design Lead — add to memory?"]
        [stale task를 갱신하고 memory gap을 채움]
```

### 업무 약어

Memory가 채워지면 Claude가 shorthand를 즉시 해석합니다.

```
You: todd한테 oracle PSR 준비해달라고 해

Claude: "Todd Martinez(Finance lead)에게 Oracle Systems deal($2.3M, Q2 close)의
         Pipeline Status Report를 준비해달라고 요청"
```

추가 질문이나 왕복 없이 처리합니다.

## 데이터 출처

> 익숙하지 않은 placeholder가 보이거나 연결된 도구를 확인해야 한다면 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

가장 좋은 경험을 위해 communication 및 project management tool을 연결하세요. 연결하지 않으면 작업과 메모리를 수동으로 관리합니다.

**포함된 MCP connection:**
- Chat(Slack): 팀 맥락 및 메시지 스캔
- Email and calendar(Microsoft 365): 액션 아이템 발견
- Knowledge base(Notion): 참고 문서
- Project tracker(Asana, Linear, Atlassian, monday.com, ClickUp): 작업 동기화
- Office suite(Microsoft 365): 문서

**추가 옵션:**
- Category별 대체 도구는 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.
