# Design 플러그인

Design productivity 플러그인입니다. 주로 Anthropic의 agentic 데스크톱 애플리케이션인 [Cowork](https://claude.com/product/cowork)용으로 설계되었지만 Claude Code에서도 동작합니다. Design critique, system management, UX writing, accessibility, research synthesis, developer handoff를 돕습니다. 어떤 design team에서도 사용할 수 있으며, 입력만으로 standalone 동작하고 Figma 및 다른 tool을 연결하면 더 강력해집니다.

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

모든 명령은 design을 설명하거나 screenshot을 붙여넣는 것만으로 **standalone** 동작하며, MCP connector를 연결하면 더 강력해집니다.

## 스킬

관련 상황에서 Claude가 자동으로 사용하는 domain knowledge입니다:

| 스킬 | 설명 |
|---|---|
| `design-critique` | Usability, hierarchy, consistency에 대한 structured design feedback을 제공합니다. "review this design", "critique this mockup", "what do you think of this screen?" 또는 exploration부터 final polish까지 어느 단계든 Figma link나 screenshot feedback 요청 시 트리거됩니다. |
| `design-system-management` | Design token, component library, pattern documentation을 관리합니다. |
| `ux-writing` | Clear, concise, consistent, brand-aligned한 effective microcopy를 작성합니다. |
| `accessibility-review` | Design 또는 page에 대해 WCAG 2.1 AA accessibility audit를 실행합니다. "audit accessibility", "check a11y", "is this accessible?" 또는 handoff 전 color contrast, keyboard navigation, touch target size, screen reader behavior를 review할 때 트리거됩니다. |
| `user-research` | User research를 계획하고 수행하며 종합합니다. "user research plan", "interview guide", "usability test", "survey design", "research questions" 또는 research를 통해 user를 이해하는 모든 단계에 도움이 필요할 때 트리거됩니다. |
| `design-handoff` | Design에서 developer handoff spec을 생성합니다. Design이 engineering 준비 상태이고 layout, design token, component prop, interaction state, responsive breakpoint, edge case, animation detail을 다루는 spec sheet가 필요할 때 사용합니다. |

## 예시 워크플로

### Design feedback 받기

```
/critique
```

Figma link, screenshot을 공유하거나 design을 설명하세요. Usability, visual hierarchy, consistency, accessibility에 대한 structured feedback을 받습니다.

### Design system audit

```
/design-system audit
```

Component library를 consistency, completeness, naming convention 기준으로 review합니다. 구체적인 improvement recommendation이 포함된 report를 받습니다.

### UX copy 작성

```
/ux-copy error messages for payment flow
```

Tone guidance, alternative, localization note가 포함된 context-appropriate copy를 받습니다.

### Developer handoff

```
/handoff
```

Figma link를 공유하면 measurement, design token, component state, interaction note, edge case가 포함된 complete spec을 받습니다.

### Accessibility check

```
/accessibility
```

Design 또는 URL을 공유하면 specific issue, severity, remediation step이 포함된 WCAG 2.1 AA compliance report를 받습니다.

### Research 종합

```
/research-synthesis
```

Interview transcript, survey result, usability test note를 upload하면 theme, insight, prioritized recommendation을 받습니다.

## 단독 사용과 통합 사용

모든 명령과 스킬은 integration 없이도 동작합니다:

| 할 수 있는 일 | 단독 사용 | 통합 사용 시 |
|-----------------|------------|-------------------|
| Design critique | 설명 또는 screenshot | Figma MCP(design 직접 수집) |
| Design system | System 설명 | Figma MCP(component library audit) |
| Handoff spec | 설명 또는 screenshot | Figma MCP(exact measurement, token) |
| UX copy | Context 설명 | Knowledge base(brand voice guideline) |
| Accessibility | 설명 또는 screenshot | Figma MCP, real usage data용 analytics |
| Research synthesis | Transcript 붙여넣기 | User feedback tool(raw data 수집) |

## MCP 통합

> 익숙하지 않은 placeholder가 보이거나 연결된 tool을 확인해야 한다면 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

더 풍부한 경험을 위해 tool을 연결하세요:

| 범주 | 예시 | 가능해지는 일 |
|---|---|---|
| **Design tool** | Figma | Design 수집, component inspect, design token 접근 |
| **User feedback** | Intercom, Productboard | Raw feedback, feature request, NPS data |
| **Project tracker** | Linear, Asana, Jira | Design을 ticket에 link하고 implementation 추적 |
| **Knowledge base** | Notion | Brand guideline, design principle, research repository |
| **Product analytics** | Amplitude, Mixpanel | Research synthesis와 design decision용 usage data |

지원되는 integration 전체 목록은 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.
