# Design 플러그인

디자인 생산성 플러그인입니다. 주로 Anthropic의 agentic 데스크톱 애플리케이션인 [Cowork](https://claude.com/product/cowork)용으로 설계되었지만 Claude Code에서도 동작합니다. 디자인 비평, 디자인 시스템 관리, UX 글쓰기, 접근성, 리서치 종합, 개발자 핸드오프를 돕습니다. 어떤 디자인 팀에서도 사용할 수 있으며, 입력만으로 단독 동작하고 Figma 및 다른 도구를 연결하면 더 강력해집니다.

## 설치

```bash
claude plugins add knowledge-work-plugins/design
```

## 명령

Slash command로 호출하는 명시적 워크플로입니다:

| 명령 | 설명 |
|---|---|
| `/critique` | 사용성, 시각적 위계, 접근성, 일관성에 대한 구조화된 디자인 피드백을 받습니다. |
| `/design-system` | 컴포넌트, 토큰, 패턴을 중심으로 디자인 시스템을 감사, 문서화, 확장합니다. |
| `/handoff` | 치수, 토큰, 상태, 상호작용, 엣지 케이스가 포함된 개발자 핸드오프 명세를 생성합니다. |
| `/ux-copy` | 마이크로카피, 오류 메시지, 빈 상태, 온보딩 플로우 등 UX 카피를 작성하거나 검토합니다. |
| `/accessibility` | WCAG 준수, 색상 대비, 스크린 리더, 키보드 탐색을 포함한 접근성 감사를 실행합니다. |
| `/research-synthesis` | 인터뷰, 설문, 사용성 테스트 등 사용자 리서치를 실행 가능한 인사이트로 종합합니다. |

모든 명령은 design을 설명하거나 screenshot을 붙여넣는 것만으로 **단독** 동작하며, MCP connector를 연결하면 더 강력해집니다.

## 스킬

관련 상황에서 Claude가 자동으로 사용하는 도메인 지식입니다:

| 스킬 | 설명 |
|---|---|
| `design-critique` | 사용성, 계층 구조, 일관성에 대한 구조화된 디자인 피드백을 제공합니다. "review this design", "critique this mockup", "what do you think of this screen?" 또는 탐색부터 최종 다듬기까지 어느 단계든 Figma 링크나 스크린샷 피드백 요청 시 트리거됩니다. |
| `design-system-management` | 디자인 토큰, 컴포넌트 라이브러리, 패턴 문서를 관리합니다. |
| `ux-writing` | 명확하고 간결하며 일관되고 브랜드에 맞는 효과적인 마이크로카피를 작성합니다. |
| `accessibility-review` | 디자인 또는 페이지에 대해 WCAG 2.1 AA 접근성 감사를 실행합니다. "audit accessibility", "check a11y", "is this accessible?" 또는 핸드오프 전 색상 대비, 키보드 탐색, 터치 대상 크기, 스크린 리더 동작을 검토할 때 트리거됩니다. |
| `user-research` | 사용자 리서치를 계획하고 수행하며 종합합니다. "user research plan", "interview guide", "usability test", "survey design", "research questions" 또는 리서치를 통해 사용자를 이해하는 모든 단계에 도움이 필요할 때 트리거됩니다. |
| `design-handoff` | 디자인에서 개발자 핸드오프 명세를 생성합니다. 디자인이 엔지니어링 준비 상태이고 레이아웃, 디자인 토큰, 컴포넌트 prop, 상호작용 상태, 반응형 breakpoint, 엣지 케이스, 애니메이션 상세를 다루는 명세서가 필요할 때 사용합니다. |

## 예시 워크플로

### Design feedback 받기

```
/critique
```

Figma link, screenshot을 공유하거나 design을 설명하세요. 사용성, 시각적 계층, 일관성, 접근성에 대한 구조화된 feedback을 받습니다.

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
