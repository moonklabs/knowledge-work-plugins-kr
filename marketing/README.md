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
| `campaign-planning` | Campaign frameworks, channel selection, content calendar creation, budget allocation, and success metrics |
| `brand-voice` | Brand voice documentation, voice attributes, tone adaptation, style guide enforcement, and terminology management |
| `competitive-analysis` | Competitive research methodology, messaging comparison, content gap analysis, positioning, and battlecard creation |
| `performance-analytics` | Key metrics by channel, reporting templates, trend analysis, attribution modeling, and optimization frameworks |

## 예시 워크플로

### Drafting a Blog Post

```
> /draft-content
Type: blog post
Topic: How AI is transforming B2B marketing
Audience: Marketing directors at mid-market SaaS companies
Key messages: AI saves time on repetitive tasks, improves personalization, requires human oversight
Tone: Authoritative but approachable
Length: 1200 words
```

Claude will generate a structured blog post draft with an engaging headline, introduction with a hook, organized sections, SEO-optimized subheadings, and a clear call to action.

### Planning a Campaign

```
> /campaign-plan
Goal: Drive 500 signups for our new product launch
Audience: Technical decision-makers at enterprise companies
Timeline: 6 weeks
Budget range: $20,000-$30,000
```

Claude will produce a campaign brief covering objectives, audience segmentation, key messages, channel strategy, a week-by-week content calendar, and KPIs to track.

### Reviewing Content Against Brand Guidelines

```
> /brand-review
[paste your draft content]
```

If your brand style guide is configured in local settings, Claude will check your content against voice, tone, terminology, and messaging pillars. If not configured, Claude will ask about your guidelines or provide a generic review for clarity, consistency, and professionalism.

## Configuration

Configure your brand voice, style guide, and target personas in a local settings file for personalized output. This allows commands like `/draft-content` and `/brand-review` to automatically apply your brand standards without prompting each time.

## MCP 통합

> If you see unfamiliar placeholders or need to check which tools are connected, see [CONNECTORS.md](CONNECTORS.md).

This plugin works with the following MCP servers:

- **Slack** — Share drafts, reports, and briefs with your team
- **Canva** — Create and edit design assets
- **Figma** — Access design files and brand assets
- **HubSpot** — Pull campaign data, manage contacts, and track marketing automation
- **Amplitude** — Pull product analytics and user behavior data for performance reporting
- **Notion** — Access briefs, style guides, and campaign documents
- **Ahrefs** — SEO keyword research, backlink analysis, and site audits
- **Similarweb** — Competitive traffic analysis and market benchmarking
- **Klaviyo** — Draft and review email marketing sequences and campaigns
- **Supermetrics** — Pull marketing data from multiple platforms for analytics and reporting
