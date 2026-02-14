# Productivity 플러그인

[Cowork](https://claude.com/product/cowork)(Anthropic의 에이전틱 데스크톱 애플리케이션)를 위해 주로 설계된 생산성 플러그인입니다. Claude Code에서도 작동합니다. 작업 관리, 업무 기억, 시각적 대시보드를 제공하며, Claude가 사람, 프로젝트, 용어를 학습하여 챗봇이 아닌 동료처럼 행동합니다.

## 설치

```
claude plugins add knowledge-work-plugins/productivity
```

## 기능

이 플러그인은 Claude에게 업무에 대한 지속적인 이해를 제공합니다:

- **작업 관리** — Claude가 읽고, 쓰고, 실행하는 마크다운 작업 목록(`TASKS.md`)입니다. 작업을 자연스럽게 추가하면 Claude가 상태를 추적하고, 오래된 항목을 분류하고, 외부 도구와 동기화합니다.
- **업무 기억** — Claude에게 약어, 사람, 프로젝트, 용어를 가르치는 2단계 기억 시스템입니다. "oracle 건에 대한 PSR을 todd에게 요청해"라고 말하면 Claude는 누가, 무엇을, 어떤 거래인지 정확히 알고 있습니다.
- **시각적 대시보드** — 작업의 보드 뷰와 Claude가 업무에 대해 알고 있는 내용의 실시간 뷰를 제공하는 로컬 HTML 파일입니다. 보드나 파일에서 편집하면 동기화가 유지됩니다.

## 커맨드

| 커맨드 | 기능 |
|--------|------|
| `/start` | 작업 및 기억 초기화, 대시보드 열기 |
| `/update` | 오래된 항목 분류, 기억 공백 확인, 외부 도구와 동기화(해당되는 경우) |
| `/update --comprehensive` | 이메일, 캘린더, 채팅 전체 스캔 — 놓친 할 일 표시 및 새로운 기억 제안 |

## 스킬

| 스킬 | 설명 |
|------|------|
| `memory-management` | 2단계 기억 시스템 — 작업 기억을 위한 CLAUDE.md, 심층 저장을 위한 memory/ 디렉토리 |
| `task-management` | 공유 TASKS.md 파일을 사용한 마크다운 기반 작업 추적 |

## 예시 워크플로우

### 시작하기

```
사용자: /start

Claude: [TASKS.md, CLAUDE.md, memory/ 디렉토리, dashboard.html 생성]
        [브라우저에서 대시보드 열기]
        [역할, 팀, 현재 우선순위에 대해 질문하여 기억 시드]
```

### 자연스럽게 작업 추가하기

```
사용자: 금요일까지 Sarah를 위한 예산 제안서를 검토하고,
       Greg와 동기화한 후 Q2 로드맵 초안을 작성하고,
       Platform 팀의 API 스펙에 대해 후속 조치를 취해야 해

Claude: [세 가지 작업 모두를 컨텍스트와 함께 TASKS.md에 추가]
        [대시보드가 자동으로 업데이트됨]
```

### 아침 동기화

```
사용자: /update --comprehensive

Claude: [이메일, 캘린더, 채팅에서 새로운 실행 항목 스캔]
        [표시: "예산 제안서 검토가 내일까지입니다 — 아직 열려 있음"]
        [제안: "3개 스레드에서 언급된 새로운 사람: Jamie Park,
         Design Lead — 기억에 추가할까요?"]
        [오래된 작업 업데이트 및 기억 공백 채우기]
```

### 업무 약어

기억이 채워지면 Claude는 약어를 즉시 해독합니다:

```
사용자: oracle 건에 대한 PSR을 todd에게 요청해

Claude: "Todd Martinez(Finance lead)에게 Oracle Systems 거래
         ($2.3M, Q2 마감)에 대한 Pipeline Status Report 준비를 요청합니다"
```

명확화 질문이나 왕복이 없습니다.

## 데이터 소스

> 익숙하지 않은 플레이스홀더가 보이거나 연결된 도구를 확인해야 하는 경우 [CONNECTORS.md](CONNECTORS.md)를 참조하세요.

최상의 경험을 위해 커뮤니케이션 및 프로젝트 관리 도구를 연결하세요. 도구가 없으면 작업과 기억을 수동으로 관리합니다.

**포함된 MCP 연결:**
- 채팅(Slack) - 팀 컨텍스트 및 메시지 스캔
- 이메일 및 캘린더(Microsoft 365) - 실행 항목 발견
- 지식 베이스(Notion) - 참조 문서
- 프로젝트 트래커(Asana, Linear, Atlassian, monday.com, ClickUp) - 작업 동기화
- Office 스위트(Microsoft 365) - 문서

**추가 옵션:**
- 각 카테고리의 대체 도구는 [CONNECTORS.md](CONNECTORS.md)를 참조하세요
