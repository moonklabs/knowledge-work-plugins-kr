# Productivity 플러그인

[Cowork](https://claude.com/product/cowork)를 위한 생산성 플러그인입니다. Anthropic의 에이전트형 데스크톱 앱인 Cowork에 맞춰 설계했지만 Claude Code에서도 동작합니다. 할 일 관리, 업무 기억, 시각 대시보드를 통해 Claude가 사람과 프로젝트, 용어를 학습해 챗봇이 아니라 동료처럼 일하도록 만듭니다.

## 설치

```
claude plugins add knowledge-work-plugins/productivity
```

## 제공 기능

이 플러그인은 Claude가 업무를 지속적으로 이해하도록 돕습니다.

- **할 일 관리**: Claude가 읽고 쓰며 실행하는 Markdown 작업 목록(`TASKS.md`)입니다. 자연스럽게 할 일을 추가하면 상태를 추적하고, 오래된 항목을 정리하며, 외부 도구와 동기화합니다.
- **업무 기억**: 단축어, 사람, 프로젝트, 용어를 Claude에게 가르치는 2단계 메모리 시스템입니다. "oracle PSR은 todd에게 부탁"이라고 말하면 누구에게 무엇을 어떤 딜에 대해 말하는지 정확히 압니다.
- **시각 대시보드**: 로컬 HTML 파일로 작업 보드와 Claude가 이해한 업무 컨텍스트를 실시간으로 보여줍니다. 보드에서 수정하든 파일에서 수정하든 동기화됩니다.

## 커맨드

| Command | 설명 |
|---------|------|
| `/start` | 작업 + 메모리를 초기화하고 대시보드를 엽니다 |
| `/update` | 오래된 항목을 정리하고 메모리의 빈틈을 확인하며 필요 시 외부 도구와 동기화합니다 |
| `/update --comprehensive` | 이메일, 캘린더, 채팅을 심층 스캔해 누락된 TODO를 표시하고 새 메모리를 제안합니다 |

## 스킬

| Skill | 설명 |
|-------|------|
| `memory-management` | 2단계 메모리 시스템: 작업 메모리용 `CLAUDE.md`, 장기 저장용 `memory/` 디렉터리 |
| `task-management` | 공유 파일 `TASKS.md`를 사용하는 Markdown 기반 작업 추적 |

## 예시 워크플로

### 시작하기

```
You: /start

Claude: [TASKS.md, CLAUDE.md, memory/ 디렉터리, dashboard.html 생성]
        [브라우저에서 대시보드 열기]
        [역할, 팀, 현재 우선순위를 물어 메모리 시드 생성]
```

### 자연스럽게 할 일 추가

```
You: 금요일까지 Sarah 예산안을 리뷰하고,
     Greg과 싱크 맞춘 뒤 Q2 로드맵을 초안으로 만들고,
     플랫폼 팀의 API 스펙에 후속 대응해야 해요

Claude: [세 가지 작업을 컨텍스트와 함께 TASKS.md에 추가]
        [대시보드 자동 업데이트]
```

### 아침 싱크

```
You: /update --comprehensive

Claude: [이메일, 캘린더, 채팅에서 새 액션 아이템 스캔]
        [표시: "예산안 리뷰는 내일까지 — 아직 미완료"]
        [제안: "3개 스레드에서 새 인물 언급: Jamie Park,
         Design Lead — 메모리에 추가할까요?"]
        [오래된 작업 정리 및 메모리 빈틈 보완]
```

### 업무 단축어

메모리가 채워지면 Claude는 단축어를 바로 해석합니다.

```
You: oracle PSR은 todd에게 부탁

Claude: "Todd Martinez(재무 리드)에게 Oracle Systems 딜에 대한
         Pipeline Status Report 작성 요청하기($2.3M, Q2 클로징)"
```

추가 질문 없이 바로 처리합니다.

## 데이터 소스

> 낯선 플레이스홀더가 보이거나 어떤 도구가 연결되어 있는지 확인하려면 [CONNECTORS.md](CONNECTORS.md)를 보세요.

최상의 경험을 위해 커뮤니케이션 및 프로젝트 관리 도구를 연결하세요. 연결하지 않아도 작업과 메모리를 수동으로 관리할 수 있습니다.

**기본 포함 MCP 연결:**
- 채팅(Slack): 팀 컨텍스트 및 메시지 스캔
- 이메일/캘린더(Microsoft 365): 액션 아이템 발견
- 지식베이스(Notion): 레퍼런스 문서
- 프로젝트 트래커(Asana, Linear, Atlassian, monday.com, ClickUp): 작업 동기화
- 오피스 제품군(Microsoft 365): 문서

**추가 옵션:**
- 카테고리별 대체 도구는 [CONNECTORS.md](CONNECTORS.md)에서 확인하세요
