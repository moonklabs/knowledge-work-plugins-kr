# 커넥터

## 도구 참조 방식

플러그인 파일은 `~~category`를 플레이스홀더로 사용하여 해당 카테고리에서 사용자가 연결한 도구를 나타냅니다. 예를 들어, `~~marketing automation`은 MCP 서버가 있는 HubSpot, Marketo 또는 기타 마케팅 플랫폼을 의미할 수 있습니다.

플러그인은 **도구에 구애받지 않으며** — 특정 제품이 아닌 카테고리(디자인, SEO, 이메일 마케팅 등) 기준으로 워크플로우를 설명합니다. `.mcp.json`에 특정 MCP 서버가 사전 구성되어 있지만, 해당 카테고리의 모든 MCP 서버가 호환됩니다.

## 이 플러그인의 커넥터

| 카테고리 | 플레이스홀더 | 포함된 서버 | 기타 옵션 |
|----------|-------------|-----------------|---------------|
| 채팅 | `~~chat` | Slack | Microsoft Teams |
| 디자인 | `~~design` | Canva, Figma | Adobe Creative Cloud |
| 마케팅 자동화 | `~~marketing automation` | HubSpot | Marketo, Pardot, Mailchimp |
| 프로덕트 분석 | `~~product analytics` | Amplitude | Mixpanel, Google Analytics |
| 지식 베이스 | `~~knowledge base` | Notion | Confluence, Guru |
| SEO | `~~SEO` | Ahrefs, Similarweb | Semrush, Moz |
| 이메일 마케팅 | `~~email marketing` | Klaviyo | Mailchimp, Brevo, Customer.io |
