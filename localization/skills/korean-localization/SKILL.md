---
name: korean-localization
description: >
  Systematically translate Claude Code plugins from English to Korean while preserving structure and technical accuracy.
  Use when users want to "한글화", "번역", "translate to Korean", "localize to Korean", "한국어로 번역", or work with Korean language content in plugins.
---

# Korean Localization

Translate Claude Code plugin files from English to Korean with precision, maintaining technical accuracy, YAML integrity, and structural consistency across all component types.

## Overview

This skill provides a systematic approach to localizing plugin content for Korean-speaking users. It handles the unique challenges of translating technical documentation while preserving:
- YAML frontmatter structure and parsing
- Code block integrity
- Technical terminology consistency
- Markdown formatting
- Plugin component relationships

**Key Principle**: Translate user-facing content; preserve technical identifiers and structure.

## Translation Scope

### What Gets Translated

| File Type | Elements to Translate | Preserve |
|-----------|----------------------|----------|
| `plugin.json` | `description` | `name`, `version`, `author` |
| `SKILL.md` | Frontmatter `description`, markdown body | Frontmatter `name`, code blocks |
| Commands (`.md`) | Frontmatter `description`, `argument-hint`, markdown body | Frontmatter `name`, `arguments` (keys) |
| README.md | All prose content | File paths, URLs, code examples |
| CONNECTORS.md | Table content, descriptions | Placeholder syntax `~~category` |
| Reference files | All prose content | Technical terms as needed |

### What Never Gets Translated

- `name` fields in all frontmatter (plugin names, skill names, command names, argument names)
- Code blocks (bash, javascript, python, json, etc.)
- File paths and directory names
- URLs and links
- Brand names (Claude, Anthropic, GitHub, Slack, etc.)
- Technical abbreviations (API, CLI, SDK, YAML, JSON, etc.)
- Variable placeholders (`${CLAUDE_PLUGIN_ROOT}`, `~~category`)
- Tool references in commands (`Read`, `Write`, `Grep`, `Bash`, etc.)

## Workflow

### Phase 1: Inventory & Planning

1. **Read the plugin.json** to understand the plugin scope
2. **List all translatable files**:
   - `plugin.json` (description field only)
   - All `SKILL.md` files in `skills/*/SKILL.md`
   - All command files in `commands/*.md`
   - `README.md`
   - `CONNECTORS.md` (if exists)
   - Reference files in `skills/*/references/*.md`
3. **Confirm scope** with the user if the file count is large (>20 files)

### Phase 2: File-by-File Translation

**For each file**, follow this sequence:

1. **Read** the original file
2. **Identify** translatable vs. non-translatable sections
3. **Translate** prose content using the glossary and guidelines
4. **Preserve** YAML structure, code blocks, and technical identifiers
5. **Validate** YAML parsing if frontmatter was modified
6. **Write** the translated file (in-place, as this is a `-kr` fork)

**Translation order** (recommended):
1. `plugin.json` (quick win, sets context)
2. `README.md` (overview understanding)
3. Skills (core content, highest value)
4. Commands (user-facing, high value)
5. Reference files (detailed content)
6. `CONNECTORS.md` (if exists)

### Phase 3: Verification & Report

After completing all translations:

1. **Check YAML integrity**: Ensure all frontmatter still parses
2. **Verify code blocks**: Confirm no code was altered
3. **Review terminology**: Check for consistency across files
4. **Report summary**:
   - Files translated
   - Word count estimates
   - Any ambiguities or decisions made

## Translation Quality Standards

### Tone & Register

- **Explanatory content**: 합쇼체 (formal polite: ~습니다/~합니다)
- **Conversational examples**: 해요체 (informal polite: ~해요/~어요)
- **Field labels & headings**: 명사형 (noun form, no verb ending)
- **Consistency**: Use the same register for similar content types across all files

### Natural Korean

- Convert English SVO to Korean SOV sentence structure
- Use Korean conjunctions (그리고, 그러나, 따라서) not literal "and", "but", "therefore"
- Prefer natural Korean phrasing over word-for-word translation
- Use 존댓말 appropriately for technical instructions

### Technical Terminology

- **Refer to `references/translation-glossary.md`** for domain-specific terms
- Keep English in parentheses on first use if helpful: "워크플로우(workflow)"
- Use established Korean tech terms: "커밋(commit)", "브랜치(branch)"
- Don't translate tool names: Slack, Linear, Jira, GitHub

### Formatting Preservation

- **Tables**: Keep column alignment, translate headers and cell content
- **Lists**: Maintain bullet/number structure, translate items
- **Headings**: Translate text, keep Markdown `#` levels
- **Links**: Translate link text, keep URLs unchanged
- **Code blocks**: Never alter code, translate preceding explanations only

## Edge Cases & Gotchas

### Mixed English/Korean in UI

When a UI element mixes technical terms with Korean:
- ✅ "GitHub 이슈(issue)를 생성합니다"
- ✅ "Linear에서 작업(task)을 추적합니다"
- ❌ "깃허브 이슈를 생성합니다" (don't translate brand names)

### YAML Frontmatter Strings

YAML strings with special characters need careful handling:

```yaml
# Original
description: Use when you want to "create a task" or "track work"

# Translated - keep quotes balanced
description: "작업 생성" 또는 "작업 추적"이 필요할 때 사용합니다
```

**Rule**: If the original has quotes, maintain them. Use YAML `>` folding for long descriptions.

### Command `argument-hint`

This field appears in the CLI autocomplete. Translate it, but keep it concise:

```yaml
# Original
argument-hint: "project name or ID"

# Translated
argument-hint: "프로젝트 이름 또는 ID"
```

### Table Cell Content

Preserve table structure; translate cell content:

```markdown
| Component | Count | Purpose |
|-----------|-------|---------|
| Skills    | 1     | Domain knowledge for X |

# Translated
| 컴포넌트 | 개수 | 목적 |
|----------|------|------|
| 스킬     | 1    | X를 위한 도메인 지식 |
```

### Code Examples in Prose

When prose references code inline, translate the prose but not the code:

```markdown
Use the `Read` tool to access files.

# Translated
파일에 접근하려면 `Read` 도구를 사용합니다.
```

### Placeholder Syntax in CONNECTORS.md

Keep `~~category` placeholders unchanged:

```markdown
Create an issue in ~~project tracker

# Translated
~~project tracker에 이슈를 생성합니다
```

## Progressive Disclosure

This SKILL.md provides the core workflow. For detailed reference material:

- **Translation Glossary**: See `references/translation-glossary.md` for 200+ domain-specific Korean terms across 11 plugin domains
- **Translation Guidelines**: See `references/translation-guidelines.md` for detailed rules on tone, sentence structure, and file-type-specific conventions

## Usage Examples

**Example 1: Translate a single plugin**

```
User: "productivity 플러그인을 한글화해줘"

Claude:
1. Reads productivity/.claude-plugin/plugin.json
2. Lists all translatable files (3 skills, 5 commands, 1 README)
3. Confirms scope with user
4. Translates files one by one
5. Reports completion with summary
```

**Example 2: Translate all plugins**

```
User: "모든 플러그인을 한국어로 번역해줘"

Claude:
1. Lists all 11 plugins from marketplace.json
2. Confirms this is ~164 files
3. Asks if user wants to proceed or start with one plugin
4. If confirmed, processes plugins sequentially
5. Reports progress and completion
```

## Quality Checklist

Before marking a file as complete:

- [ ] YAML frontmatter parses correctly (if applicable)
- [ ] All code blocks are unchanged
- [ ] `name` fields are unchanged
- [ ] Table structure is preserved
- [ ] Technical terms match glossary
- [ ] Tone is consistent with file type
- [ ] No English remains in translatable sections (except technical terms)

## Additional Resources

- `references/translation-glossary.md` — 11 domain-specific terminology tables
- `references/translation-guidelines.md` — Detailed translation rules and conventions
