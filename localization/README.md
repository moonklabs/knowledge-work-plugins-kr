# Localization Plugin

> Systematically translate Claude Code plugins from English to Korean while preserving structure, technical accuracy, and YAML integrity.

> Claude Code 플러그인을 영문에서 한국어로 체계적으로 번역하며, 구조, 기술 정확성, YAML 무결성을 유지합니다.

## Overview | 개요

The Localization plugin provides a structured workflow for translating all components of Claude Code plugins into Korean. It handles the unique challenges of technical translation including:

- YAML frontmatter preservation
- Code block integrity
- Technical terminology consistency
- Markdown structure maintenance
- Domain-specific vocabulary

이 플러그인은 Claude Code 플러그인의 모든 구성 요소를 한국어로 번역하는 구조화된 워크플로우를 제공합니다. 기술 번역의 고유한 과제를 처리합니다:

- YAML 프론트매터 보존
- 코드 블록 무결성
- 기술 용어 일관성
- 마크다운 구조 유지
- 도메인별 어휘

## Components | 구성 요소

### Skills

- **korean-localization**: Core translation workflow with phase-by-phase guidance
  - 단계별 가이드를 제공하는 핵심 번역 워크플로우

### Reference Files

Located in `skills/korean-localization/references/`:

- **translation-glossary.md**: 200+ domain-specific Korean terms across 11 plugin domains
  - 11개 플러그인 도메인에 걸친 200개 이상의 도메인별 한국어 용어

- **translation-guidelines.md**: Detailed rules for tone, sentence structure, and file-type-specific conventions
  - 어조, 문장 구조, 파일 유형별 규칙에 대한 상세 가이드

## Usage | 사용법

### Basic Translation | 기본 번역

To translate a plugin, use natural language triggers:

```
"productivity 플러그인을 한글화해줘"
"Translate the sales plugin to Korean"
"legal 플러그인을 번역해줘"
```

### Translation Workflow | 번역 워크플로우

The skill follows a three-phase workflow:

1. **Inventory & Planning**: Lists all translatable files and confirms scope
   - **인벤토리 및 계획**: 번역 가능한 모든 파일을 나열하고 범위를 확인

2. **File-by-File Translation**: Translates each file while preserving structure
   - **파일별 번역**: 구조를 유지하면서 각 파일을 번역

3. **Verification & Report**: Validates YAML, checks code integrity, reports completion
   - **검증 및 리포트**: YAML 검증, 코드 무결성 확인, 완료 보고

### What Gets Translated | 번역 대상

- ✅ `plugin.json` description field
- ✅ SKILL.md frontmatter descriptions
- ✅ SKILL.md markdown body
- ✅ Command frontmatter (description, argument-hint)
- ✅ Command markdown body
- ✅ README.md prose content
- ✅ CONNECTORS.md content
- ✅ Reference files

### What Never Gets Translated | 번역 제외 항목

- ❌ `name` fields (plugins, skills, commands, arguments)
- ❌ Code blocks (bash, python, javascript, etc.)
- ❌ Brand names (Claude, Slack, GitHub, etc.)
- ❌ Tool names (Read, Write, Grep, Bash)
- ❌ File paths and URLs
- ❌ Technical abbreviations (API, CLI, SDK, YAML)

## Translation Quality Standards | 번역 품질 기준

### Tone & Register | 어조 및 격식

- **합쇼체** (formal polite): Technical explanations, documentation
- **해요체** (informal polite): Conversational examples, dialogue
- **명사형** (noun form): Headings, labels, field names

### Terminology Consistency | 용어 일관성

All technical terms follow the glossary to ensure consistency across:
- 11 plugin domains (Customer Support, Sales, Finance, Legal, Marketing, Product Management, Data, Enterprise Search, Bio Research, Productivity, and Common terms)
- 200+ standardized translations
- Clear guidelines on when to keep English vs. translate

모든 기술 용어는 다음 전반에 걸쳐 일관성을 보장하기 위해 용어집을 따릅니다:
- 11개 플러그인 도메인
- 200개 이상의 표준화된 번역
- 영문 유지 vs. 번역에 대한 명확한 가이드라인

## Examples | 예시

### Example 1: Single Plugin Translation

```
User: "customer-support 플러그인을 한글화해줘"

Claude:
1. Reads customer-support/.claude-plugin/plugin.json
2. Lists all translatable files (skills, commands, README)
3. Confirms scope with user
4. Translates files sequentially
5. Validates YAML and code integrity
6. Reports completion summary
```

### Example 2: Specific File Translation

```
User: "sales/skills/cold-outreach/SKILL.md를 번역해줘"

Claude:
1. Reads the SKILL.md file
2. Translates frontmatter description
3. Translates markdown body
4. Preserves code blocks and technical identifiers
5. Validates YAML parsing
6. Writes translated file
```

## Domain Coverage | 도메인 범위

The glossary covers terminology from these domains:

1. **Common Terms**: plugin, skill, command, workflow, task
2. **Customer Support**: ticket, triage, escalation, knowledge base
3. **Sales**: prospect, pipeline, outreach, deal
4. **Finance**: journal entry, reconciliation, audit
5. **Legal**: contract review, NDA, compliance
6. **Marketing**: campaign, SEO, conversion, engagement
7. **Product Management**: feature spec, roadmap, user story
8. **Data**: SQL query, dataset, visualization
9. **Enterprise Search**: knowledge synthesis, semantic search
10. **Bio Research**: genomics, RNA-seq, clinical trial
11. **Productivity**: task management, calendar, reminder

## File Structure | 파일 구조

```
localization/
├── .claude-plugin/
│   └── plugin.json              # Plugin manifest
├── README.md                    # This file
└── skills/
    └── korean-localization/
        ├── SKILL.md             # Core translation workflow
        └── references/
            ├── translation-glossary.md      # Domain-specific terms
            └── translation-guidelines.md    # Translation rules
```

## Installation | 설치

This plugin is part of the `knowledge-work-plugins-kr` repository. To use it:

이 플러그인은 `knowledge-work-plugins-kr` 레포지토리의 일부입니다. 사용하려면:

```bash
# If using from the forked repository
cd knowledge-work-plugins-kr/localization

# The plugin is automatically available in the repository context
# 레포지토리 컨텍스트에서 자동으로 사용 가능합니다
```

## Contributing | 기여

When adding new domain terms to the glossary:

1. Add to the appropriate domain section in `translation-glossary.md`
2. Include English term, Korean translation, and usage notes
3. Maintain alphabetical order within each section
4. Update this README if adding a new domain section

## License

Same as the parent `knowledge-work-plugins` repository.

---

**Version**: 0.1.0
**Author**: Anthropic
**Repository**: knowledge-work-plugins-kr
