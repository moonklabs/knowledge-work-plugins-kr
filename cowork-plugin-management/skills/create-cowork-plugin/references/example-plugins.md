# 플러그인 예시

다양한 복잡도 수준의 완전한 플러그인 구조 3가지입니다. 4단계에서 구현할 때 템플릿으로 사용합니다.

## 최소 플러그인: 단일 커맨드

슬래시 커맨드 하나만 있고 다른 컴포넌트가 없는 간단한 플러그인입니다.

### 구조

```
meeting-notes/
├── .claude-plugin/
│   └── plugin.json
├── commands/
│   └── meeting-notes.md
└── README.md
```

### plugin.json

```json
{
  "name": "meeting-notes",
  "version": "0.1.0",
  "description": "전사 기록(transcript)에서 구조화된 회의록 생성",
  "author": {
    "name": "사용자"
  }
}
```

### commands/meeting-notes.md

```markdown
---
description: 전사 기록에서 구조화된 회의록 생성
argument-hint: [전사-기록-파일]
allowed-tools: Read, Write
---

@$1에 있는 전사 기록을 읽고 구조화된 회의록을 생성하세요.

다음 섹션을 포함해야 합니다:
1. **참석자** — 언급된 모든 참가자 목록
2. **요약** — 회의에 대한 2-3문장 중심의 개요
3. **주요 의사결정** — 내려진 결정 사항에 대한 번호 매기기 목록
4. **실행 항목** — 소유자, 작업, 마감일 컬럼이 포함된 테이블
5. **미결 질문** — 해결되지 않은 모든 사안

회의록을 전사 기록 파일명 뒤에 `-notes`를 붙인 새 파일로 저장하세요.
```

---

## 표준 플러그인: 스킬 + 커맨드 + MCP

도메인 지식, 사용자 커맨드, 외부 서비스 통합을 결합한 플러그인입니다.

### 구조

```
code-quality/
├── .claude-plugin/
│   └── plugin.json
├── commands/
│   ├── review.md
│   └── fix-lint.md
├── skills/
│   └── coding-standards/
│       ├── SKILL.md
│       └── references/
│           └── style-rules.md
├── .mcp.json
└── README.md
```

### plugin.json

```json
{
  "name": "code-quality",
  "version": "0.1.0",
  "description": "리뷰, 린팅 및 스타일 가이드를 통한 코딩 표준 준수",
  "author": {
    "name": "사용자"
  }
}
```

### commands/review.md

```markdown
---
description: 코드 변경 사항의 스타일 및 품질 이슈 리뷰
allowed-tools: Read, Grep, Bash(git:*)
---

변경 된 파일 목록 가져오기: !`git diff --name-only`

각 변경된 파일에 대해:
1. 파일을 읽습니다
2. coding-standards 스킬을 참고하여 스타일 위반 사항을 확인합니다
3. 잠재적 버그 또는 안티 패턴을 식별합니다
4. 보안 우려 사항을 플래그 지정합니다

다음을 포함하는 요약 보고서를 제시합니다:
- 파일 경로
- 이슈 심각도 (Error, Warning, Info)
- 설명 및 수정 제안
```

### commands/fix-lint.md

```markdown
---
description: 변경된 파일의 린팅 이슈 자동 수정
allowed-tools: Read, Write, Edit, Bash(npm:*)
---

린터 실행: !`npm run lint -- --format json 2>&1`

린터 출력을 파싱하고 각 이슈를 수정합니다:
- 자동 수정이 가능한 이슈는 수정을 직접 적용합니다
- 수동 수정이 필요한 이슈는 프로젝트 컨벤션에 따라 수정합니다
- 아키텍처 변경이 필요한 이슈는 건너뜁니다

모든 수정이 끝나면 린터를 다시 실행하여 깨끗한 출력을 확인합니다.
```

### skills/coding-standards/SKILL.md

```yaml
---
name: coding-standards
description: >
  사용자가 "코딩 표준", "스타일 가이드", "명명 규칙", "코드 포맷팅 규칙"에 대해 묻거나
  프로젝트 고유의 코드 품질 기대치에 대한 가이드가 필요할 때 이 스킬을 사용합니다.
version: 0.1.0
---
```

```markdown
# 코딩 표준 (Coding Standards)

일관되고 고품질인 코드를 위한 프로젝트 코딩 표준 및 컨벤션입니다.

## 핵심 규칙

- 변수와 함수에는 camelCase를 사용합니다
- 클래스와 타입에는 PascalCase를 사용합니다
- let 대신 const를 선호하고, var는 사용하지 않습니다
- 최대 라인 길이: 100자
- 모든 노출(exported) 함수에는 명시적 반환 타입을 사용합니다

## 임포트 순서

1. 외부 패키지
2. 내부 패키지 (@/ 에일리어스 사용)
3. 상대 경로 임포트
4. 타입 전용 임포트는 마지막에

## 추가 리소스

- **`references/style-rules.md`** — 언어별 전체 스타일 규칙
```

### .mcp.json

```json
{
  "mcpServers": {
    "github": {
      "type": "http",
      "url": "https://api.githubcopilot.com/mcp/"
    }
  }
}
```

---

## 전체 기능 플러그인: 모든 컴포넌트 유형

스킬, 커맨드, 에이전트, 훅, MCP 통합을 도구 비종속적 커넥터와 함께 사용하는 플러그인입니다.

### 구조

```
engineering-workflow/
├── .claude-plugin/
│   └── plugin.json
├── commands/
│   ├── standup-prep.md
│   └── create-ticket.md
├── skills/
│   └── team-processes/
│       ├── SKILL.md
│       └── references/
│           └── workflow-guide.md
├── agents/
│   └── ticket-analyzer.md
├── hooks/
│   └── hooks.json
├── .mcp.json
├── CONNECTORS.md
└── README.md
```

### plugin.json

```json
{
  "name": "engineering-workflow",
  "version": "0.1.0",
  "description": "엔지니어링 워크플로우 간소화: 스탠드업 준비, 티켓 관리 및 코드 품질",
  "author": {
    "name": "사용자"
  },
  "keywords": ["engineering", "workflow", "tickets", "standup"]
}
```

### agents/ticket-analyzer.md

```markdown
---
name: ticket-analyzer
description: 사용자가 티켓을 분석하거나, 들어오는 이슈를 분류하거나, 백로그의 우선순위를 정해야 할 때 이 에이전트를 사용합니다.

<example>
컨텍스트: 사용자가 스프린트 계획을 준비 중임
사용자: "이 새 티켓들의 분류를 도와줘"
어시스턴트: "티켓을 검토하고 분류하기 위해 ticket-analyzer 에이전트를 사용하겠습니다."
<commentary>
티켓 분류는 여러 차원에 걸친 체계적인 분석이 필요하므로 에이전트가 적합합니다.
</commentary>
</example>

<example>
컨텍스트: 사용자가 방대한 백로그를 가지고 있음
사용자: "다음 스프린트를 위해 내 백로그의 우선순위를 정해줘"
어시스턴트: "백로그를 분석하여 우선순위를 제안하기 위해 ticket-analyzer 에이전트를 가동하겠습니다."
<commentary>
백로그 우선순위 지정은 에이전트에 적합한 다단계 자율 작업입니다.
</commentary>
</example>

model: inherit
color: cyan
tools: ["Read", "Grep"]
---

당신은 티켓 분석 전문가입니다. 티켓의 우선순위, 작업량 및 의존성을 분석하십시오.

**핵심 책임:**
1. 티켓을 유형별로 분류 (버그, 기능, 기술 부채, 개선)
2. 상대적 작업량 추정 (S, M, L, XL)
3. 티켓 간의 의존성 식별
4. 우선순위 순서 권장

**분석 프로세스:**
1. 모든 티켓 설명 읽기
2. 각 티켓을 유형별로 분류
3. 범위에 따른 작업량 추정
4. 의존성 맵 작성
5. 영향도 대비 작업량 비율로 순위 지정

**출력 형식:**
| 티켓 | 유형 | 작업량 | 의존성 | 우선순위 |
|------|------|--------|--------|----------|
| ...  | ...  | ...    | ...    | ...      |

이어서 상위 5개 우선순위에 대한 간략한 근거를 제시합니다.
```

### hooks/hooks.json

```json
{
  "SessionStart": [
    {
      "matcher": "",
      "hooks": [
        {
          "type": "command",
          "command": "echo '## 팀 컨텍스트\n\n스프린트 주기: 2주. 스탠업: 매일 오전 9:30. 티켓 관리는 ~~project tracker를 사용하세요.'",
          "timeout": 5
        }
      ]
    }
  ]
}
```

### CONNECTORS.md

```markdown
# 커넥터 (Connectors)

## 도구 참조 작동 방식

플러그인 파일은 사용자가 해당 카테고리에 연결한 도구의 플레이스홀더로 `~~카테고리`를 사용합니다. 플러그인은 도구에 구애받지 않습니다.

## 이 플러그인의 커넥터

| 카테고리 | 플레이스홀더 | 포함된 서버 | 기타 옵션 |
|----------|-------------|-----------------|---------------|
| 프로젝트 트래커 | `~~project tracker` | Linear | Asana, Jira, Monday |
| 채팅 | `~~chat` | Slack | Microsoft Teams |
| 소스 제어 | `~~source control` | GitHub | GitLab, Bitbucket |
```

### .mcp.json

```json
{
  "mcpServers": {
    "linear": {
      "type": "sse",
      "url": "https://mcp.linear.app/sse"
    },
    "github": {
      "type": "http",
      "url": "https://api.githubcopilot.com/mcp/"
    },
    "slack": {
      "type": "http",
      "url": "https://slack.mcp.claude.com/mcp"
    }
  }
}
```
