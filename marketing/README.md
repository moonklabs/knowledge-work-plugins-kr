# Marketing 플러그인

Anthropic의 에이전트 데스크톱 애플리케이션인 [Cowork](https://claude.com/product/cowork)을 위해 설계된 마케팅 플러그인으로, Claude Code에서도 사용할 수 있습니다. 콘텐츠 제작, 캠페인 기획, 브랜드 보이스 관리, 경쟁사 분석, 성과 리포팅 기능을 제공합니다.

## 설치

```bash
claude plugins add knowledge-work-plugins/marketing
```

## 커맨드

| 커맨드 | 설명 |
|---|---|
| `/draft-content` | 블로그 포스트, 소셜 미디어, 이메일 뉴스레터, 랜딩 페이지, 보도자료, 사례 연구를 작성합니다 |
| `/campaign-plan` | 목표, 채널, 콘텐츠 캘린더, 성공 지표를 포함한 전체 캠페인 브리프를 생성합니다 |
| `/brand-review` | 브랜드 보이스, 스타일 가이드, 메시징 핵심 요소에 맞춰 콘텐츠를 검토합니다 |
| `/competitive-brief` | 경쟁사를 조사하고 포지셔닝 및 메시징 비교 자료를 생성합니다 |
| `/performance-report` | 핵심 지표, 트렌드, 최적화 권장 사항을 포함한 마케팅 성과 리포트를 작성합니다 |
| `/seo-audit` | 키워드 리서치, 온페이지 분석, 콘텐츠 갭, 기술 점검, 경쟁사 비교를 포함한 종합 SEO 감사를 수행합니다 |
| `/email-sequence` | 너처 플로우, 온보딩, 드립 캠페인 등을 위한 다중 이메일 시퀀스를 설계하고 작성합니다 |

## 스킬

| 스킬 | 설명 |
|---|---|
| `content-creation` | 콘텐츠 유형별 템플릿, 채널별 작성 모범 사례, SEO 기초, 헤드라인 공식, CTA 가이드를 제공합니다 |
| `campaign-planning` | 캠페인 프레임워크, 채널 선택, 콘텐츠 캘린더 생성, 예산 배분, 성공 지표를 다룹니다 |
| `brand-voice` | 브랜드 보이스 문서화, 보이스 속성, 톤 조정, 스타일 가이드 적용, 용어 관리를 지원합니다 |
| `competitive-analysis` | 경쟁사 리서치 방법론, 메시징 비교, 콘텐츠 갭 분석, 포지셔닝, 배틀카드 작성을 제공합니다 |
| `performance-analytics` | 채널별 핵심 지표, 리포팅 템플릿, 트렌드 분석, 기여도 모델링, 최적화 프레임워크를 다룹니다 |

## 워크플로우 예시

### 블로그 포스트 작성

```
> /draft-content
Type: blog post
Topic: How AI is transforming B2B marketing
Audience: Marketing directors at mid-market SaaS companies
Key messages: AI saves time on repetitive tasks, improves personalization, requires human oversight
Tone: Authoritative but approachable
Length: 1200 words
```

Claude가 매력적인 헤드라인, 후크가 포함된 서론, 체계적으로 정리된 섹션, SEO에 최적화된 소제목, 명확한 CTA(Call to Action)를 갖춘 구조화된 블로그 포스트 초안을 생성합니다.

### 캠페인 기획

```
> /campaign-plan
Goal: Drive 500 signups for our new product launch
Audience: Technical decision-makers at enterprise companies
Timeline: 6 weeks
Budget range: $20,000-$30,000
```

Claude가 목표, 오디언스 세분화, 핵심 메시지, 채널 전략, 주간 콘텐츠 캘린더, 추적할 KPI를 포함한 캠페인 브리프를 작성합니다.

### 브랜드 가이드라인에 맞춘 콘텐츠 검토

```
> /brand-review
[paste your draft content]
```

브랜드 스타일 가이드가 로컬 설정에 구성되어 있으면, Claude가 보이스, 톤, 용어, 메시징 핵심 요소에 맞춰 콘텐츠를 점검합니다. 구성되어 있지 않은 경우, Claude가 가이드라인에 대해 질문하거나 명확성, 일관성, 전문성 관점에서 일반적인 검토를 제공합니다.

## 구성

로컬 설정 파일에서 브랜드 보이스, 스타일 가이드, 타겟 페르소나를 구성하면 개인화된 결과물을 얻을 수 있습니다. 이를 통해 `/draft-content`, `/brand-review` 등의 커맨드가 매번 입력하지 않아도 자동으로 브랜드 기준을 적용합니다.

## MCP 통합

> 익숙하지 않은 플레이스홀더가 보이거나 연결된 도구를 확인하려면 [CONNECTORS.md](CONNECTORS.md)를 참조하세요.

이 플러그인은 다음 MCP 서버와 연동됩니다:

- **Slack** — 팀에 초안, 리포트, 브리프를 공유합니다
- **Canva** — 디자인 에셋을 생성하고 편집합니다
- **Figma** — 디자인 파일과 브랜드 에셋에 접근합니다
- **HubSpot** — 캠페인 데이터를 조회하고, 연락처를 관리하며, 마케팅 자동화를 추적합니다
- **Amplitude** — 성과 리포팅을 위한 프로덕트 분석 및 사용자 행동 데이터를 조회합니다
- **Notion** — 브리프, 스타일 가이드, 캠페인 문서에 접근합니다
- **Ahrefs** — SEO 키워드 리서치, 백링크 분석, 사이트 감사를 수행합니다
- **Similarweb** — 경쟁사 트래픽 분석과 시장 벤치마킹을 수행합니다
- **Klaviyo** — 이메일 마케팅 시퀀스와 캠페인을 작성하고 검토합니다
