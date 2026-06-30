# Marketing 플러그인

Marketing 플러그인입니다. 주로 Anthropic의 agentic 데스크톱 애플리케이션인 [Cowork](https://claude.com/product/cowork)용으로 설계되었지만 Claude Code에서도 동작합니다. 콘텐츠 제작, 캠페인 계획, 브랜드 보이스 관리, 경쟁 분석, 성과 보고를 지원합니다.

## 설치

```bash
claude plugins add knowledge-work-plugins/marketing
```

## 명령

| 명령 | 설명 |
|---|---|
| `/draft-content` | 블로그 글, 소셜 게시물, 이메일 뉴스레터, 랜딩 페이지, 보도자료, 사례 연구 초안을 작성합니다. |
| `/campaign-plan` | 목표, 채널, 콘텐츠 캘린더, 성공 지표가 포함된 전체 캠페인 브리프를 생성합니다. |
| `/brand-review` | 브랜드 보이스, 스타일 가이드, 메시징 원칙에 맞춰 콘텐츠를 검토합니다. |
| `/competitive-brief` | 경쟁사를 조사하고 포지셔닝과 메시지 비교를 생성합니다. |
| `/performance-report` | 핵심 지표, 추세, 최적화 추천이 포함된 마케팅 성과 보고서를 만듭니다. |
| `/seo-audit` | 키워드 리서치, 온페이지 분석, 콘텐츠 갭, 기술 점검, 경쟁사 비교를 포함한 SEO 감사를 실행합니다. |
| `/email-sequence` | 리드 육성, 온보딩, 드립 캠페인을 위한 다중 이메일 시퀀스를 설계하고 초안을 작성합니다. |

## 스킬

| 스킬 | 설명 |
|---|---|
| `content-creation` | 여러 채널의 마케팅 콘텐츠 초안을 작성합니다. 채널별 서식, SEO 최적화 문구, headline option, call to action이 필요할 때 사용합니다. |
| `campaign-planning` | 캠페인 프레임워크, 채널 선택, 콘텐츠 캘린더, 예산 배분, 성공 지표를 다룹니다. |
| `brand-voice` | 브랜드 보이스 문서화, voice attribute, tone adaptation, style guide enforcement, 용어 관리를 다룹니다. |
| `competitive-analysis` | 경쟁 리서치 방법론, 메시지 비교, 콘텐츠 갭 분석, 포지셔닝, battlecard 작성을 다룹니다. |
| `performance-analytics` | 채널별 핵심 지표, 보고 템플릿, 추세 분석, attribution modeling, 최적화 프레임워크를 다룹니다. |

## 예시 워크플로

### Blog post 초안 작성

```text
> /draft-content
Type: blog post
Topic: How AI is transforming B2B marketing
Audience: Marketing directors at mid-market SaaS companies
Key messages: AI saves time on repetitive tasks, improves personalization, requires human oversight
Tone: Authoritative but approachable
Length: 1200 words
```

Claude는 눈길을 끄는 headline, hook이 있는 introduction, 정리된 section, SEO-optimized subheading, 명확한 call to action이 포함된 structured blog post draft를 생성합니다.

### Campaign 계획

```text
> /campaign-plan
Goal: 새 product launch에서 signup 500건 달성
Audience: Technical decision-makers at enterprise companies
Timeline: 6 weeks
Budget range: $20,000-$30,000
```

Claude는 objective, audience segmentation, key message, channel strategy, week-by-week content calendar, tracking KPI를 다루는 campaign brief를 생성합니다.

### Brand guideline 기준 content review

```text
> /brand-review
[paste your draft content]
```

Local settings에 brand style guide가 설정되어 있으면 Claude가 voice, tone, terminology, messaging pillar 기준으로 content를 확인합니다. 설정되어 있지 않으면 guideline을 묻거나 clarity, consistency, professionalism 기준의 generic review를 제공합니다.

## 설정

Personalized output을 위해 local settings file에 brand voice, style guide, target persona를 설정하세요. 이렇게 하면 `/draft-content`, `/brand-review` 같은 command가 매번 묻지 않고 brand standard를 자동 적용합니다.

## MCP 통합

> 익숙하지 않은 placeholder가 보이거나 연결된 도구를 확인해야 한다면 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

이 플러그인은 다음 MCP server와 함께 동작합니다.

- **Slack**: Draft, report, brief를 team과 공유
- **Canva**: Design asset 생성 및 편집
- **Figma**: Design file 및 brand asset 접근
- **HubSpot**: Campaign data 수집, contact 관리, marketing automation tracking
- **Amplitude**: Performance reporting용 product analytics 및 user behavior data 수집
- **Notion**: Brief, style guide, campaign document 접근
- **Ahrefs**: SEO keyword research, backlink analysis, site audit
- **Similarweb**: Competitive traffic analysis 및 market benchmarking
- **Klaviyo**: Email marketing sequence와 campaign 초안 작성 및 review
- **Supermetrics**: Analytics/reporting용 multi-platform marketing data 수집
