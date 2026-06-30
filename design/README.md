# Design 플러그인

Design productivity plugin입니다. 주로 Anthropic의 agentic desktop application인 [Cowork](https://claude.com/product/cowork)용으로 설계되었지만 Claude Code에서도 동작합니다. Design critique, system management, UX writing, accessibility, research synthesis, developer handoff를 돕습니다. 어떤 design team에서도 사용할 수 있으며, 입력만으로 standalone 동작하고 Figma 및 다른 tool을 연결하면 더 강력해집니다.

## 설치

```bash
claude plugins add knowledge-work-plugins/design
```

## 명령

Slash command로 호출하는 명시적 workflow입니다:

| 명령 | 설명 |
|---|---|
| `/critique` | Usability, visual hierarchy, accessibility, consistency에 대한 structured design feedback을 받습니다. |
| `/design-system` | Component, token, pattern을 중심으로 design system을 audit, document, extend합니다. |
| `/handoff` | Measurement, token, state, interaction, edge case가 포함된 developer handoff spec을 생성합니다. |
| `/ux-copy` | Microcopy, error message, empty state, onboarding flow 등 UX copy를 작성하거나 review합니다. |
| `/accessibility` | WCAG compliance, color contrast, screen reader, keyboard navigation을 포함한 accessibility audit를 실행합니다. |
| `/research-synthesis` | Interview, survey, usability test 등 user research를 actionable insight로 종합합니다. |

All commands work **standalone** (describe your design or paste screenshots) and get **supercharged** with MCP connectors.

## 스킬

관련 상황에서 Claude가 자동으로 사용하는 domain knowledge입니다:

| 스킬 | 설명 |
|---|---|
| `design-critique` | Usability, hierarchy, consistency에 대한 structured design feedback을 제공합니다. "review this design", "critique this mockup", "what do you think of this screen?" 또는 exploration부터 final polish까지 어느 단계든 Figma link나 screenshot feedback 요청 시 트리거됩니다. |
| `design-system-management` | Manage design tokens, component libraries, and pattern documentation |
| `ux-writing` | Write effective microcopy — clear, concise, consistent, and brand-aligned |
| `accessibility-review` | Design 또는 page에 대해 WCAG 2.1 AA accessibility audit를 실행합니다. "audit accessibility", "check a11y", "is this accessible?" 또는 handoff 전 color contrast, keyboard navigation, touch target size, screen reader behavior를 review할 때 트리거됩니다. |
| `user-research` | User research를 계획하고 수행하며 종합합니다. "user research plan", "interview guide", "usability test", "survey design", "research questions" 또는 research를 통해 user를 이해하는 모든 단계에 도움이 필요할 때 트리거됩니다. |
| `design-handoff` | Design에서 developer handoff spec을 생성합니다. Design이 engineering 준비 상태이고 layout, design token, component prop, interaction state, responsive breakpoint, edge case, animation detail을 다루는 spec sheet가 필요할 때 사용합니다. |

## 예시 워크플로

### Getting Design Feedback

```
/critique
```

Share a Figma link, screenshot, or describe your design. Get structured feedback on usability, visual hierarchy, consistency, and accessibility.

### Auditing Your Design System

```
/design-system audit
```

I'll review your component library for consistency, completeness, and naming conventions. Get a report with specific improvement recommendations.

### Writing UX Copy

```
/ux-copy error messages for payment flow
```

Get context-appropriate copy with tone guidance, alternatives, and localization notes.

### Developer Handoff

```
/handoff
```

Share a Figma link and get a complete spec: measurements, design tokens, component states, interaction notes, and edge cases.

### Accessibility Check

```
/accessibility
```

Share a design or URL. Get a WCAG 2.1 AA compliance report with specific issues, severity, and remediation steps.

### Synthesizing Research

```
/research-synthesis
```

Upload interview transcripts, survey results, or usability test notes. Get themes, insights, and prioritized recommendations.

## 단독 사용과 통합 사용

모든 명령과 스킬은 integration 없이도 동작합니다:

| 할 수 있는 일 | 단독 사용 | 통합 사용 시 |
|-----------------|------------|-------------------|
| Design critique | Describe or screenshot | Figma MCP (pull designs directly) |
| Design system | Describe your system | Figma MCP (audit component library) |
| Handoff specs | Describe or screenshot | Figma MCP (exact measurements, tokens) |
| UX copy | Describe the context | Knowledge base (brand voice guidelines) |
| Accessibility | Describe or screenshot | Figma MCP, analytics for real usage data |
| Research synthesis | Paste transcripts | User feedback tools (pull raw data) |

## MCP 통합

> If you see unfamiliar placeholders or need to check which tools are connected, see [CONNECTORS.md](CONNECTORS.md).

더 풍부한 경험을 위해 tool을 연결하세요:

| 범주 | 예시 | 가능해지는 일 |
|---|---|---|
| **Design tool** | Figma | Pull designs, inspect components, access design tokens |
| **User feedback** | Intercom, Productboard | Raw feedback, feature requests, NPS data |
| **Project tracker** | Linear, Asana, Jira | Link designs to tickets, track implementation |
| **Knowledge base** | Notion | Brand guidelines, design principles, research repository |
| **Product analytics** | Amplitude, Mixpanel | Usage data for research synthesis and design decisions |

지원되는 integration 전체 목록은 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.
