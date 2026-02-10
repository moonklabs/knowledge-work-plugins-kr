# Translation Guidelines

Detailed rules for translating Claude Code plugin files from English to Korean while preserving technical accuracy, structure, and YAML integrity.

## Table of Contents

1. [Tone & Register](#tone--register)
2. [Sentence Structure](#sentence-structure)
3. [File-Type-Specific Rules](#file-type-specific-rules)
4. [Prohibited Translations](#prohibited-translations)
5. [Quality Checklist](#quality-checklist)

---

## Tone & Register

Korean has multiple levels of formality (존댓말). Use the appropriate register for each content type.

### 합쇼체 (Formal Polite: ~습니다/~합니다)

**Use for**: Technical explanations, instructions, documentation

**Characteristics**:
- Most formal polite form
- Professional and authoritative
- Used in formal writing, technical docs, manuals

**Examples**:
```markdown
# English
This skill helps you manage tasks efficiently.

# Korean (합쇼체)
이 스킬은 작업을 효율적으로 관리하는 데 도움을 줍니다.
```

```markdown
# English
Use the Read tool to access files.

# Korean (합쇼체)
파일에 접근하려면 Read 도구를 사용합니다.
```

**When to use**: README.md, SKILL.md body, reference files, technical explanations

### 해요체 (Informal Polite: ~해요/~어요)

**Use for**: Conversational examples, dialogue, user-facing messages

**Characteristics**:
- Polite but friendly
- Conversational tone
- Used in casual professional communication

**Examples**:
```markdown
# English (in a conversational example)
User: "Can you help me with this?"
Claude: "Sure! Let me check the documentation."

# Korean (해요체)
사용자: "이거 도와줄 수 있어요?"
Claude: "네! 문서를 확인해 볼게요."
```

**When to use**: Dialogue examples, user quotes, conversational contexts

### 명사형 (Noun Form)

**Use for**: Headings, labels, field names, table headers

**Characteristics**:
- No verb ending
- Concise and clear
- Used for navigation and labels

**Examples**:
```markdown
# English
## Getting Started

# Korean (명사형)
## 시작하기
```

```yaml
# English
description: Task management system

# Korean (명사형)
description: 작업 관리 시스템
```

**When to use**: Headings, plugin.json fields, table headers, navigation items

### Register Consistency

**Rule**: Use the same register for similar content types across all files.

| Content Type | Register | Example |
|--------------|----------|---------|
| Plugin descriptions | 명사형 | "작업 관리 시스템" |
| Skill body | 합쇼체 | "이 스킬은 작업을 관리합니다" |
| README sections | 합쇼체 | "플러그인을 설치합니다" |
| Dialogue examples | 해요체 | "이거 도와줄게요" |
| Table headers | 명사형 | "설명", "목적" |
| Bullet points | 합쇼체 | "작업을 생성합니다" |

---

## Sentence Structure

English and Korean have fundamentally different word orders. Translate meaning, not word-for-word.

### SVO → SOV Conversion

**English**: Subject-Verb-Object (SVO)
**Korean**: Subject-Object-Verb (SOV)

```markdown
# English (SVO)
This plugin manages your tasks efficiently.

# Korean (SOV)
이 플러그인은 작업을 효율적으로 관리합니다.
      ↑       ↑          ↑
   Subject  Object    Verb
```

### Conjunctions

Use natural Korean conjunctions, not literal translations.

| English | ❌ Literal | ✅ Natural Korean |
|---------|-----------|------------------|
| and | 그리고 | 및, ~와/과, ~고 |
| but | 하지만 | 그러나, ~지만 |
| so, therefore | 그래서 | 따라서, ~므로 |
| if | 만약 | ~면, ~라면 |
| when | 언제 | ~때, ~할 때 |
| because | 왜냐하면 | ~기 때문에, ~므로 |

**Examples**:

```markdown
# English
If you want to create a task, use this command.

# ❌ Literal
만약 당신이 작업을 생성하고 싶으면, 이 커맨드를 사용하세요.

# ✅ Natural
작업을 생성하려면 이 커맨드를 사용합니다.
```

```markdown
# English
This skill is powerful but requires setup.

# ❌ Literal
이 스킬은 강력하다 하지만 설정이 필요하다.

# ✅ Natural
이 스킬은 강력하지만 설정이 필요합니다.
```

### Omitting Subjects

Korean often omits subjects when context is clear. Don't force "당신(you)" into translations.

```markdown
# English
You can use this skill to manage tasks.

# ❌ Too literal
당신은 작업을 관리하기 위해 이 스킬을 사용할 수 있습니다.

# ✅ Natural (subject omitted)
이 스킬을 사용하여 작업을 관리할 수 있습니다.
```

### Compound Sentences

Break complex English sentences into natural Korean clauses.

```markdown
# English
This command searches your codebase for files matching a pattern, filters them by type, and returns the top results.

# ❌ Single run-on sentence
이 커맨드는 패턴과 일치하는 파일을 코드베이스에서 검색하고, 타입별로 필터링하고, 상위 결과를 반환합니다.

# ✅ Natural with better flow
이 커맨드는 패턴과 일치하는 파일을 코드베이스에서 검색합니다. 파일을 타입별로 필터링한 후 상위 결과를 반환합니다.
```

---

## File-Type-Specific Rules

Different file types have different requirements. Follow these conventions.

### plugin.json

**Translate**: `description` field only
**Preserve**: `name`, `version`, `author`, all keys

```json
{
  "name": "productivity",  // ❌ DO NOT TRANSLATE
  "version": "0.1.0",      // ❌ DO NOT TRANSLATE
  "description": "작업 관리, 일정 계획, 중요한 업무 맥락 기억을 지원합니다. 캘린더, 이메일, 채팅과 동기화하여 모든 것을 정리하고 추적합니다.",  // ✅ TRANSLATE
  "author": {
    "name": "Anthropic"    // ❌ DO NOT TRANSLATE
  }
}
```

**Register**: 명사형 or 합쇼체 (no strict rule, match existing style)

### SKILL.md Frontmatter

**YAML structure**: Must remain valid YAML
**Translate**: `description` field
**Preserve**: `name` field

```yaml
---
name: task-management  # ❌ DO NOT TRANSLATE
description: >         # ✅ TRANSLATE (use > for multi-line)
  작업을 관리하고, 우선순위를 지정하고, 마감일을 추적합니다.
  "작업 관리", "할 일 추가", "작업 추적"이 필요할 때 사용합니다.
---
```

**Special handling for quotes**:

```yaml
# English
description: Use when you want to "create a task" or "track work"

# Korean - maintain quote balance
description: "작업 생성" 또는 "작업 추적"이 필요할 때 사용합니다
```

**Multi-line with `>`**:

```yaml
# Use > for folding long descriptions
description: >
  이 스킬은 작업을 관리하고, 일정을 계획하고,
  중요한 업무 맥락을 기억하는 데 도움을 줍니다.
  캘린더, 이메일, 채팅과 동기화하여 모든 것을 정리합니다.
```

### Command Frontmatter

**Translate**: `description`, `argument-hint`
**Preserve**: `name`, `arguments` (keys), `tools`

```yaml
---
name: create-task       # ❌ DO NOT TRANSLATE
description: 새로운 작업을 생성하고 할 일 목록에 추가합니다  # ✅ TRANSLATE
argument-hint: "작업 제목 및 설명"  # ✅ TRANSLATE (keep concise)
arguments:
  - name: title         # ❌ DO NOT TRANSLATE
    description: 작업 제목  # ✅ TRANSLATE
  - name: due_date      # ❌ DO NOT TRANSLATE
    description: 마감일 (YYYY-MM-DD)  # ✅ TRANSLATE
tools:                  # ❌ DO NOT TRANSLATE
  - Read
  - Write
---
```

**Register**: 합쇼체 for descriptions, 명사형 for argument-hint

### Markdown Body

**Translate**: All prose content
**Preserve**: Code blocks, file paths, URLs, technical identifiers

```markdown
# English
## Overview

This skill helps you manage tasks. Use the `Read` tool to access task files.

Available tools:
- `Read`
- `Write`
- `Grep`

Example usage:
```bash
claude task create "Review PR"
```

# Korean
## 개요

이 스킬은 작업 관리를 돕습니다. 작업 파일에 접근하려면 `Read` 도구를 사용합니다.

사용 가능한 도구:
- `Read`
- `Write`
- `Grep`

사용 예시:
```bash
claude task create "Review PR"
```
```

**Key points**:
- Translate headings (## 개요)
- Translate prose sentences
- Keep code blocks unchanged (```bash ... ```)
- Keep inline code unchanged (`Read`)
- Translate list items that are prose, preserve technical terms

### Code Blocks

**Rule**: NEVER translate code blocks

```markdown
# ✅ CORRECT
```python
def create_task(title):
    return {"title": title, "status": "pending"}
```

# ❌ INCORRECT (DO NOT DO THIS)
```python
def 작업_생성(제목):
    return {"제목": 제목, "상태": "대기중"}
```
```

**Translate surrounding explanations**, not code:

```markdown
# English
Use this Python function to create a task:
```python
def create_task(title):
    return {"title": title}
```

# Korean
작업을 생성하려면 다음 Python 함수를 사용합니다:
```python
def create_task(title):
    return {"title": title}
```
```

### Tables

**Translate**: Headers and cell content
**Preserve**: Structure (column alignment, `|` separators)

```markdown
# English
| Component | Count | Purpose |
|-----------|-------|---------|
| Skills    | 3     | Domain knowledge |
| Commands  | 5     | User actions |

# Korean
| 컴포넌트 | 개수 | 목적 |
|----------|------|------|
| 스킬     | 3    | 도메인 지식 |
| 커맨드   | 5    | 사용자 액션 |
```

**Alignment**: Keep column widths reasonable for Korean characters (wider than English)

### CONNECTORS.md

**Translate**: Prose, table content, descriptions
**Preserve**: Placeholder syntax (`~~category`)

```markdown
# English
Plugin files use `~~category` as a placeholder for whatever tool the user connects.

| Category | Placeholder | Options |
|----------|-------------|---------|
| Chat | `~~chat` | Slack, Microsoft Teams |

# Korean
플러그인 파일은 사용자가 연결한 도구를 나타내기 위해 `~~category`를 플레이스홀더로 사용합니다.

| 카테고리 | 플레이스홀더 | 옵션 |
|----------|-------------|------|
| 채팅     | `~~chat`    | Slack, Microsoft Teams |
```

**Key**: `~~chat`, `~~project tracker` etc. are NEVER translated

### README.md

**Translate**: All prose
**Preserve**: Code examples, file paths, URLs, installation commands

```markdown
# English
# Productivity Plugin

Manage tasks, plan your day, and build memory of important work context.

## Installation

```bash
claude plugin install productivity
```

## Usage

Use `/task` to create a new task.

# Korean
# Productivity 플러그인

작업 관리, 일정 계획, 중요한 업무 맥락 기억을 지원합니다.

## 설치

```bash
claude plugin install productivity
```

## 사용법

새로운 작업을 생성하려면 `/task`를 사용합니다.
```

**Register**: 합쇼체 for body, 명사형 for headings

---

## Prohibited Translations

These must NEVER be translated.

### Brand Names

Keep all brand names in English:
- ✅ Slack, Linear, Jira, GitHub, Notion, Asana
- ✅ Claude, Anthropic, OpenAI
- ✅ Google, Microsoft, Apple

### Technical Abbreviations

Keep widely-used tech abbreviations in English:
- ✅ API, CLI, SDK, REST, GraphQL
- ✅ JSON, YAML, XML, HTML, CSS
- ✅ HTTP, HTTPS, SSH, SSL, TLS
- ✅ JWT, OAuth, SAML
- ✅ CRUD (Create, Read, Update, Delete)
- ✅ SQL (in database context - but "SQL" the sales term can be translated in sales plugin)

### Claude Tool Names

Keep Claude's tool names in English (these are technical identifiers):
- ✅ Read, Write, Edit, Grep, Glob, Bash
- ✅ Task, WebFetch, WebSearch

These appear as inline code in commands and should never be translated.

### File Paths and Variables

Never translate:
- File paths: `./commands/*.md`, `/path/to/file`
- Environment variables: `${CLAUDE_PLUGIN_ROOT}`, `$HOME`
- Placeholders: `~~project tracker`, `~~chat`
- URLs: `https://example.com`
- Git refs: `main`, `master`, `HEAD`

### Plugin Component Names

Never translate the `name` field in any frontmatter:

```yaml
# ❌ WRONG
name: 작업-관리
description: 작업을 관리합니다

# ✅ CORRECT
name: task-management
description: 작업을 관리합니다
```

This applies to:
- Plugin names in `plugin.json`
- Skill names in `SKILL.md`
- Command names in command files
- Argument names in command frontmatter
- Agent names (if used)

---

## Quality Checklist

Before marking a file as complete, verify:

### YAML Integrity

- [ ] Frontmatter still parses as valid YAML (test with a YAML parser if unsure)
- [ ] Quotes are balanced (opening `"` has closing `"`)
- [ ] Multi-line strings use `>` or `|` correctly
- [ ] No unterminated strings

### Code Preservation

- [ ] All code blocks are unchanged (```language ... ```)
- [ ] Inline code is unchanged (`tool_name`)
- [ ] No Korean characters inside code blocks or backticks

### Name Fields

- [ ] All `name` fields in frontmatter are unchanged
- [ ] Plugin names, skill names, command names are in English
- [ ] Argument names are in English

### Structure

- [ ] Table structure is preserved (same number of columns, `|` alignment)
- [ ] Heading levels are unchanged (`#`, `##`, `###`)
- [ ] List formatting is preserved (bullet `-` or number `1.`)
- [ ] File paths are unchanged

### Terminology

- [ ] Technical terms match the glossary
- [ ] Same English term is translated the same way throughout the file
- [ ] Brand names are in English
- [ ] Tool names are in English

### Tone Consistency

- [ ] Register is appropriate for content type
- [ ] Consistent register throughout the file
- [ ] No awkward literal translations (check for 당신, 만약, etc.)

### Completeness

- [ ] All translatable prose is translated
- [ ] No English remains in sentences (except technical terms)
- [ ] Parenthetical English is included for technical terms on first use

---

## Common Mistakes to Avoid

### 1. Translating `name` Fields

```yaml
# ❌ WRONG
name: 작업-생성
description: 작업을 생성합니다

# ✅ CORRECT
name: create-task
description: 작업을 생성합니다
```

### 2. Altering Code Blocks

```markdown
# ❌ WRONG
작업을 생성하려면:
```python
def 작업_생성(제목):
    return {"제목": 제목}
```

# ✅ CORRECT
작업을 생성하려면:
```python
def create_task(title):
    return {"title": title}
```
```

### 3. Breaking YAML Strings

```yaml
# ❌ WRONG (missing closing quote)
description: "작업을 생성하고 관리합니다

# ✅ CORRECT
description: "작업을 생성하고 관리합니다"
```

### 4. Inconsistent Terminology

```markdown
# ❌ WRONG (inconsistent)
작업을 생성합니다. 그런 다음 태스크를 관리합니다.

# ✅ CORRECT (consistent)
작업을 생성합니다. 그런 다음 작업을 관리합니다.
```

### 5. Too Literal Translation

```markdown
# English
If you want to create a task, use this command.

# ❌ Too literal
만약 당신이 작업을 생성하기를 원한다면, 이 커맨드를 사용하세요.

# ✅ Natural
작업을 생성하려면 이 커맨드를 사용합니다.
```

### 6. Translating Tool Names

```markdown
# ❌ WRONG
파일을 읽으려면 읽기 도구를 사용합니다.

# ✅ CORRECT
파일을 읽으려면 Read 도구를 사용합니다.
```

### 7. Translating Placeholders

```markdown
# ❌ WRONG
~~프로젝트_트래커에 이슈를 생성합니다

# ✅ CORRECT
~~project tracker에 이슈를 생성합니다
```

---

## Final Notes

### When in Doubt

1. **Check the glossary first** for standard translations
2. **Preserve structure** over literal translation
3. **Prioritize naturalness** in Korean over word-for-word mapping
4. **Test YAML parsing** if you modified frontmatter
5. **Be consistent** with previous translations in the same file

### Progressive Improvement

First pass: Focus on accuracy and structure preservation
Second pass: Refine naturalness and tone
Third pass: Check consistency and terminology

### Getting Help

If you encounter:
- Ambiguous technical terms → Check glossary, ask user if unclear
- Complex YAML → Validate with a YAML parser
- Domain-specific jargon → Consult domain section of glossary
- Unclear context → Read surrounding files for clues
