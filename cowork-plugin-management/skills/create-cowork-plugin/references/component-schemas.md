# 컴포넌트 스키마

모든 플러그인 컴포넌트 유형에 대한 상세 형식 명세입니다. 4단계에서 컴포넌트를 구현할 때 참조합니다.

## 커맨드

**위치**: `commands/command-name.md`
**형식**: 선택적 YAML 프론트매터가 포함된 Markdown

### 프론트매터 필드

| 필드 | 필수 | 타입 | 설명 |
|------|------|------|------|
| `description` | 아니오 | String | `/help`에 표시되는 간단한 설명 (60자 미만) |
| `allowed-tools` | 아니오 | String 또는 Array | 커맨드가 사용할 수 있는 도구 |
| `model` | 아니오 | String | 모델 재지정: `sonnet`, `opus`, `haiku` |
| `argument-hint` | 아니오 | String | 자동완성을 위한 예상 인수 설명 |

### 커맨드 예시

```markdown
---
description: 코드의 보안 이슈 리뷰
allowed-tools: Read, Grep, Bash(git:*)
argument-hint: [파일-경로]
---

다음 항목을 포함하여 @$1의 보안 취약점을 리뷰해 주세요:
- SQL 인젝션
- XSS 공격
- 인증 우회
- 안전하지 않은 데이터 처리

구체적인 라인 번호, 심각도 등급 및 수정 제안을 제공해 주세요.
```

### 주요 규칙

- **커맨드는 Claude에 대한 지시사항입니다**, 사용자를 위한 메시지가 아닙니다. 지시문으로 작성합니다.
- `$ARGUMENTS`는 모든 인수를 단일 문자열로 캡처하고, `$1`, `$2`, `$3`은 위치별 인수를 캡처합니다.
- `@path` 구문은 커맨드 컨텍스트에 파일 내용을 포함합니다.
- `!` 백틱 구문은 동적 컨텍스트를 위해 bash를 인라인으로 실행합니다(예: `` !`git diff --name-only` ``).
- 플러그인 파일을 이식 가능하게 참조하려면 `${CLAUDE_PLUGIN_ROOT}`를 사용합니다.

### allowed-tools 패턴

```yaml
# Specific tools
allowed-tools: Read, Write, Edit, Bash(git:*)

# Bash with specific commands only
allowed-tools: Bash(npm:*), Read

# MCP tools (specific)
allowed-tools: ["mcp__plugin_name_server__tool_name"]
```

## 스킬

**위치**: `skills/skill-name/SKILL.md`
**형식**: YAML 프론트매터가 포함된 Markdown

### 프론트매터 필드

| 필드 | 필수 | 타입 | 설명 |
|------|------|------|------|
| `name` | 예 | String | 스킬 식별자 |
| `description` | 예 | String | 트리거 문구가 포함된 3인칭 설명 |
| `version` | 아니오 | String | 시맨틱 버전 |

### 스킬 예시

```yaml
---
name: api-design
description: >
  사용자가 "API 설계", "API 엔드포인트 생성", "API 구조 리뷰"를 요청하거나
  REST API 모범 사례, 엔드포인트 명명 규칙, 요청/응답 설계에 대한 가이드가
  필요할 때 이 스킬을 사용합니다.
version: 0.1.0
---
```

### 작성 스타일 규칙

- **프론트매터 설명**: 3인칭("This skill should be used when..."), 따옴표 안에 구체적인 트리거 문구를 포함합니다.
- **본문**: 명령형/부정사 형태("설정 파일을 파싱합니다", "설정 파일을 파싱해야 합니다"가 아님).
- **길이**: SKILL.md 본문을 3,000단어 미만으로 유지합니다(이상적으로는 1,500-2,000단어). 상세 내용은 `references/`로 이동합니다.

### 스킬 디렉토리 구조

```
skill-name/
├── SKILL.md              # Core knowledge (required)
├── references/           # Detailed docs loaded on demand
│   ├── patterns.md
│   └── advanced.md
├── examples/             # Working code examples
│   └── sample-config.json
└── scripts/              # Utility scripts
    └── validate.sh
```

### 점진적 공개 수준

1. **메타데이터** (항상 컨텍스트에 포함): 이름 + 설명 (~100단어)
2. **SKILL.md 본문** (스킬 트리거 시): 핵심 지식 (<5k 단어)
3. **번들 리소스** (필요 시): 참조, 예시, 스크립트 (무제한)

## 에이전트

**위치**: `agents/agent-name.md`
**형식**: YAML 프론트매터가 포함된 Markdown

### 프론트매터 필드

| 필드 | 필수 | 타입 | 설명 |
|------|------|------|------|
| `name` | 예 | String | 소문자, 하이픈, 3-50자 |
| `description` | 예 | String | `<example>` 블록이 포함된 트리거 조건 |
| `model` | 예 | String | `inherit`, `sonnet`, `opus`, 또는 `haiku` |
| `color` | 예 | String | `blue`, `cyan`, `green`, `yellow`, `magenta`, `red` |
| `tools` | 아니오 | Array | 특정 도구로 제한 |

### 에이전트 예시

```markdown
---
name: code-reviewer
description: 사용자가 철저한 코드 리뷰를 요청하거나 코드 품질, 보안 및 모범 사례에 대한 자세한 분석을 원할 때 이 에이전트를 사용합니다.

<example>
컨텍스트: 사용자가 방금 새 모듈을 작성함
사용자: "이 코드에 대해 심층 리뷰를 해 줄 수 있어?"
어시스턴트: "철저한 분석을 위해 code-reviewer 에이전트를 사용하겠습니다."
<commentary>
사용자가 이 에이전트의 전문 분야인 상세 리뷰를 명시적으로 요청했습니다.
</commentary>
</example>

<example>
컨텍스트: 사용자가 PR을 머지하려 함
사용자: "머지하기 전에 이거 리뷰해 줘"
어시스턴트: "code-reviewer 에이전트를 사용하여 종합적인 리뷰를 실행하겠습니다."
<commentary>
머지 전 리뷰는 에이전트의 구조화된 분석 프로세스를 통해 도움을 받을 수 있습니다.
</commentary>
</example>

model: inherit
color: blue
tools: ["Read", "Grep", "Glob"]
---

당신은 보안, 성능, 유지보수성 및 정확성 전반에 걸쳐 이슈를 식별하는 데 집중하는 코드 리뷰 전문가입니다.

**핵심 책임:**
1. 코드 구조 및 조직 분석
2. 보안 취약점 식별
3. 성능 관련 우려 사항 플래그 지정
4. 모범 사례 준수 여부 확인

**분석 프로세스:**
1. 범위 내의 모든 파일 읽기
2. 패턴 및 안티 패턴 식별
3. 발견 사항을 심각도별로 분류
4. 구체적인 수정 제안 제공

**출력 형식:**
발견 사항을 심각도별(Critical, Warning, Info)로 그룹화하여 다음 내용을 포함하여 제시합니다:
- 파일 경로 및 라인 번호
- 이슈 설명
- 제안하는 수정 방법
```

### 에이전트 네이밍 규칙

- 3-50자
- 소문자, 숫자, 하이픈만 사용
- 영숫자로 시작하고 끝나야 함
- 밑줄, 공백, 특수 문자 사용 불가

### 색상 가이드라인

- Blue/Cyan: 분석, 리뷰
- Green: 성공 지향 작업
- Yellow: 주의, 검증
- Red: 중요, 보안
- Magenta: 창의적, 생성

## 훅

**위치**: `hooks/hooks.json`
**형식**: JSON

### 사용 가능한 이벤트

| 이벤트 | 실행 시점 |
|--------|-----------|
| `PreToolUse` | 도구 호출이 실행되기 전 |
| `PostToolUse` | 도구 호출이 완료된 후 |
| `Stop` | Claude가 응답을 완료했을 때 |
| `SubagentStop` | 하위 에이전트가 완료했을 때 |
| `SessionStart` | 세션이 시작될 때 |
| `SessionEnd` | 세션이 종료될 때 |
| `UserPromptSubmit` | 사용자가 메시지를 보낼 때 |
| `PreCompact` | 컨텍스트 압축 전 |
| `Notification` | 알림이 발생할 때 |

### 훅 유형

**프롬프트 기반** (복잡한 로직에 권장):
```json
{
  "type": "prompt",
  "prompt": "이 파일 쓰기가 프로젝트 컨벤션을 따르는지 평가하세요: $TOOL_INPUT",
  "timeout": 30
}
```
지원되는 이벤트: Stop, SubagentStop, UserPromptSubmit, PreToolUse.

**커맨드 기반** (결정론적 검사):
```json
{
  "type": "command",
  "command": "bash ${CLAUDE_PLUGIN_ROOT}/hooks/scripts/validate.sh",
  "timeout": 60
}
```

### hooks.json 예시

```json
{
  "PreToolUse": [
    {
      "matcher": "Write|Edit",
      "hooks": [
        {
          "type": "prompt",
          "prompt": "이 파일 쓰기가 프로젝트 코딩 표준을 준수하는지 확인하십시오. 표준을 위반하는 경우 이유를 설명하고 차단하십시오.",
          "timeout": 30
        }
      ]
    }
  ],
  "SessionStart": [
    {
      "matcher": "",
      "hooks": [
        {
          "type": "command",
          "command": "cat ${CLAUDE_PLUGIN_ROOT}/context/project-context.md",
          "timeout": 10
        }
      ]
    }
  ]
}
```

### 훅 출력 형식 (커맨드 훅)

커맨드 훅은 stdout으로 JSON을 반환합니다:

```json
{
  "decision": "block",
  "reason": "File write violates naming convention"
}
```

결정 값: `approve`, `block`, `ask_user` (확인 요청).

## MCP 서버

**위치**: 플러그인 루트의 `.mcp.json`
**형식**: JSON

### 서버 유형

**stdio** (로컬 프로세스):
```json
{
  "mcpServers": {
    "my-server": {
      "command": "node",
      "args": ["${CLAUDE_PLUGIN_ROOT}/servers/server.js"],
      "env": {
        "API_KEY": "${API_KEY}"
      }
    }
  }
}
```

**SSE** (OAuth를 사용하는 호스팅):
```json
{
  "mcpServers": {
    "asana": {
      "type": "sse",
      "url": "https://mcp.asana.com/sse"
    }
  }
}
```

**HTTP** (REST API):
```json
{
  "mcpServers": {
    "api-service": {
      "type": "http",
      "url": "https://api.example.com/mcp",
      "headers": {
        "Authorization": "Bearer ${API_TOKEN}"
      }
    }
  }
}
```

### 서버 유형 선택 가이드

| 유형 | 최적 용도 | 인증 방법 |
|------|-----------|-----------|
| stdio | 로컬 도구, 사용자 지정 서버 | 환경 변수 |
| SSE | 호스팅 클라우드 서비스 | OAuth (자동) |
| HTTP | REST API 백엔드 | 토큰 헤더 |

### 환경 변수 확장

모든 MCP 설정은 `${VAR_NAME}` 치환을 지원합니다:
- `${CLAUDE_PLUGIN_ROOT}` — 플러그인 디렉토리 (이식성을 위해 항상 사용)
- `${ANY_ENV_VAR}` — 사용자 환경 변수

모든 필수 환경 변수를 플러그인 README에 문서화합니다.

## CONNECTORS.md

**위치**: 플러그인 루트
**생성 시기**: 플러그인이 특정 제품이 아닌 카테고리로 외부 도구를 참조하는 경우

### 형식

```markdown
# 커넥터 (Connectors)

## 도구 참조 작동 방식

플러그인 파일은 사용자가 해당 카테고리에 연결한 도구의 플레이스홀더로 `~~카테고리`를 사용합니다. 예를 들어, `~~project tracker`는 Asana, Linear, Jira 또는 MCP 서버가 있는 기타 프로젝트 트래커를 의미할 수 있습니다.

플러그인은 도구에 구애받지 않습니다 — 즉, 특정 제품이 아닌 카테고리 관점에서 워크플로우를 서술합니다.

## 이 플러그인의 커넥터

| 카테고리 | 플레이스홀더 | 포함된 서버 | 기타 옵션 |
|----------|-------------|-----------------|---------------|
| 채팅 | `~~chat` | Slack | Microsoft Teams, Discord |
| 프로젝트 트래커 | `~~project tracker` | Linear | Asana, Jira, Monday |
```

### ~~ 플레이스홀더 사용

플러그인 파일(커맨드, 스킬, 에이전트)에서 도구를 일반적으로 참조합니다:

```markdown
Check ~~project tracker for open tickets assigned to the user.
Post a summary to ~~chat in the team channel.
```

커스터마이징 시(cowork-plugin-customizer 스킬을 통해) 이러한 참조가 특정 도구 이름으로 교체됩니다.

## README.md

모든 플러그인에는 다음을 포함하는 README가 있어야 합니다:

1. **개요** — 플러그인이 하는 일
2. **컴포넌트** — 커맨드, 스킬, 에이전트, 훅, MCP 서버 목록
3. **설정** — 필수 환경 변수 또는 설정
4. **사용법** — 각 커맨드 사용 방법 또는 각 스킬 트리거 방법
5. **커스터마이징** — CONNECTORS.md가 있으면 언급
