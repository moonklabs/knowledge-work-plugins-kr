# Marketing 플러그인

Marketing plugin입니다. 주로 Anthropic의 agentic desktop application인 [Cowork](https://claude.com/product/cowork)용으로 설계되었지만 Claude Code에서도 동작합니다. Content creation, campaign planning, brand voice management, competitive analysis, performance reporting을 지원합니다.

## 설치

```bash
claude plugins add knowledge-work-plugins/marketing
```

## 명령

| 명령 | 설명 |
|---|---|
| `/draft-content` | Blog post, social media, email newsletter, landing page, press release, case study 초안을 작성합니다. |
| `/campaign-plan` | Objective, channel, content calendar, success metric이 포함된 full campaign brief를 생성합니다. |
| `/brand-review` | Brand voice, style guide, messaging pillar 기준으로 content를 review합니다. |
| `/competitive-brief` | Competitor를 research하고 positioning 및 messaging comparison을 생성합니다. |
| `/performance-report` | Key metric, trend, optimization recommendation이 포함된 marketing performance report를 만듭니다. |
| `/seo-audit` | Keyword research, on-page analysis, content gap, technical check, competitor comparison을 포함한 SEO audit를 실행합니다. |
| `/email-sequence` | Nurture flow, onboarding, drip campaign 등을 위한 multi-email sequence를 설계하고 초안을 작성합니다. |

## 스킬

| 스킬 | 설명 |
|---|---|
| `content-creation` | Blog post, social media, email newsletter, landing page, press release, case study 등 channel 전반의 marketing content 초안을 작성합니다. Marketing content 작성, channel-specific formatting, SEO-optimized copy, headline option, call to action이 필요할 때 사용합니다. |
| `campaign-planning` | Campaign framework, channel selection, content calendar creation, budget allocation, success metric을 다룹니다. |
| `brand-voice` | Brand voice documentation, voice attribute, tone adaptation, style guide enforcement, terminology management를 다룹니다. |
| `competitive-analysis` | Competitive research methodology, messaging comparison, content gap analysis, positioning, battlecard creation을 다룹니다. |
| `performance-analytics` | Channel별 key metric, reporting template, trend analysis, attribution modeling, optimization framework를 다룹니다. |

## 예시 워크플로

### Blog post 초안 작성

```
> /draft-content
Type: blog post
Topic: How AI is transforming B2B marketing
Audience: Marketing directors at mid-market SaaS companies
Key messages: AI saves time on repetitive tasks, improves personalization, requires human oversight
Tone: Authoritative but approachable
Length: 1200 words
```

Claude는 engaging headline, hook이 있는 introduction, organized section, SEO-optimized subheading, clear call to action이 포함된 structured blog post draft를 생성합니다.

### Campaign 계획

```
> /campaign-plan
Goal: Drive 500 signups for our new product launch
Audience: Technical decision-makers at enterprise companies
Timeline: 6 weeks
Budget range: $20,000-$30,000
```

Claude는 objective, audience segmentation, key message, channel strategy, week-by-week content calendar, tracking KPI를 다루는 campaign brief를 생성합니다.

### Brand guideline 기준 content review

```
> /brand-review
[paste your draft content]
```

Local settings에 brand style guide가 설정되어 있으면 Claude가 voice, tone, terminology, messaging pillar 기준으로 content를 확인합니다. 설정되어 있지 않으면 guideline을 묻거나 clarity, consistency, professionalism 기준의 generic review를 제공합니다.

## 설정

Personalized output을 위해 local settings file에 brand voice, style guide, target persona를 설정하세요. 이렇게 하면 `/draft-content`, `/brand-review` 같은 command가 매번 묻지 않고 brand standard를 자동 적용합니다.

## MCP 통합

> 익숙하지 않은 placeholder가 보이거나 연결된 tool을 확인해야 한다면 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

이 plugin은 다음 MCP server와 함께 동작합니다.

- **Slack** — Draft, report, brief를 team과 공유
- **Canva** — Design asset 생성 및 편집
- **Figma** — Design file 및 brand asset 접근
- **HubSpot** — Campaign data 수집, contact 관리, marketing automation tracking
- **Amplitude** — Performance reporting용 product analytics 및 user behavior data 수집
- **Notion** — Brief, style guide, campaign document 접근
- **Ahrefs** — SEO keyword research, backlink analysis, site audit
- **Similarweb** — Competitive traffic analysis 및 market benchmarking
- **Klaviyo** — Email marketing sequence와 campaign 초안 작성 및 review
- **Supermetrics** — Analytics/reporting용 multi-platform marketing data 수집
