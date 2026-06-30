# Legal Productivity 플러그인

사내 법무팀을 위한 AI 기반 생산성 플러그인입니다. 주로 Anthropic의 agentic 데스크톱 애플리케이션인 [Cowork](https://claude.com/product/cowork)용으로 설계되었지만 Claude Code에서도 동작합니다. 계약 검토, NDA 분류, 컴플라이언스 워크플로, 법무 브리핑, 템플릿 응답을 자동화하며, 조직의 업무 기준과 위험 허용도에 맞게 모두 설정할 수 있습니다.

> **면책 고지:** 이 플러그인은 법무 워크플로를 보조하지만 법률 자문을 제공하지 않습니다. 결론은 항상 자격을 갖춘 법률 전문가와 검증하세요. AI가 생성한 분석은 법적 의사결정에 사용하기 전에 면허 있는 변호사의 검토를 받아야 합니다. 이 플러그인의 기본 playbook 예시는 미국 법적 입장과 관할(Delaware, New York, California)을 반영합니다. 다른 법체계(EU, UK, Netherlands, Australia 등)에서 운영한다면, 플러그인 분석에 의존하기 전에 `.claude/legal.local.md`의 playbook을 해당 관할의 법적 요구사항, 표준 계약 조건, 컴플라이언스 의무에 맞게 커스터마이즈해야 합니다.

## 대상 사용자

- **Commercial Counsel** -- 계약 협상, 벤더 관리, 거래 지원을 담당하는 법무 담당자
- **Product Counsel** -- 제품 검토, 이용약관, 개인정보 처리방침, IP 사안을 다루는 법무 담당자
- **Privacy / Compliance** -- 데이터 보호 규제, DPA 검토, 정보주체 요청, 규제 모니터링 담당자
- **Litigation Support** -- 증거보존 통지, 문서 검토 준비, 사건 브리핑 지원 담당자

## 설치

```
claude plugins add knowledge-work-plugins/legal
```

## 빠른 시작

### 1. 플러그인 설치

```
claude plugins add knowledge-work-plugins/legal
```

### 2. Playbook 설정

조직의 표준 입장을 정의하는 로컬 설정 파일을 만드세요. 팀의 협상 기준, 위험 허용도, 표준 계약 조건을 여기에 기록합니다.

Claude가 찾을 수 있는 위치에 `legal.local.md` 파일을 만드세요.

- **Cowork**: Folder picker로 Cowork와 공유한 아무 폴더에 저장합니다. 플러그인이 자동으로 찾습니다.
- **Claude Code**: 프로젝트의 `.claude/` 디렉터리에 저장합니다.

```markdown
# Legal Playbook 설정

## 계약 검토 입장

### 책임 제한
- 표준 입장: 양방향 책임 한도를 지급/지급 예정 수수료 12개월분으로 제한
- 허용 범위: 수수료 6~24개월분
- Escalation trigger: 무제한 책임, consequential damages 포함

### 면책
- 표준 입장: IP 침해와 데이터 침해에 대한 상호 면책
- 허용 가능: 제3자 청구로 제한된 면책
- Escalation trigger: 일방적 면책 의무, 무제한 면책

### IP 소유권
- 표준 입장: 각 당사자는 기존 IP를 보유하고 customer data는 customer가 소유
- Escalation trigger: 광범위한 IP 양도 조항, 기존 IP에 대한 work-for-hire 조항

### 데이터 보호
- 표준 입장: 모든 개인정보 처리에 DPA 요구
- 요구사항: 하위 처리자 통지, 종료 시 데이터 삭제, 72시간 내 침해 통지
- Escalation trigger: DPA 미제공, 보호조치 없는 국외 이전

### 기간 및 해지
- 표준 입장: 연 단위 계약, 30일 전 통지로 임의 해지 가능
- 허용 가능: 최초 기간 이후 임의 해지가 가능한 다년 계약
- Escalation trigger: 통지 기간 없는 자동 갱신, 임의 해지 불가

### 준거법
- 선호: [Your jurisdiction]
- 허용 가능: 주요 상사 관할(NY, DE, CA, England & Wales)
- Escalation trigger: 비표준 관할, 불리한 장소의 강제 중재

## NDA 기본값
- 상호 의무 필요
- 기간: 표준 2~3년, 영업비밀은 5년
- 표준 예외: 독자 개발, 공개 정보, 제3자로부터 정당하게 수령한 정보
- Residuals 조항: 범위가 좁으면 허용 가능

## 응답 템플릿
일반 문의용 template file 경로를 설정하거나 inline template을 정의하세요.
```

### 3. 도구 연결

이 플러그인은 MCP를 통해 기존 도구에 연결할 때 가장 잘 동작합니다. 미리 구성된 서버에는 Slack, Box, Egnyte, Atlassian, Microsoft 365가 포함됩니다. 지원 범주와 선택지 전체 목록은 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

## 명령

### `/review-contract` -- Playbook 기반 계약 검토

조직의 협상 기준으로 계약을 검토합니다. 표준과 다른 조항을 표시하고 redline을 생성하며 사업 영향 분석을 제공합니다.

```
/review-contract
```

입력: 파일 업로드, URL, 붙여 넣은 계약 텍스트. 맥락(your side, deadline, focus areas)을 물은 뒤 설정된 playbook 기준으로 조항별 검토를 수행합니다.

### `/triage-nda` -- NDA 사전 분류

수신 NDA를 표준 기준으로 빠르게 분류합니다. GREEN(표준 승인), YELLOW(법무 검토), RED(중요 이슈)로 나눕니다.

```
/triage-nda
```

### `/vendor-check` -- 벤더 계약 상태

연결된 시스템 전반에서 특정 벤더와의 기존 계약 상태를 확인합니다.

```
/vendor-check [vendor name]
```

기존 NDA, MSA, DPA, 만료일, 핵심 조건을 보고합니다.

### `/brief` -- 법무팀 브리핑

법무 업무를 위한 맥락 브리핑을 생성합니다.

```
/brief daily          # 법무 관련 항목의 아침 브리프
/brief topic [query]  # 특정 법무 질문에 대한 리서치 브리프
/brief incident       # 진행 중인 상황에 대한 빠른 브리프
```

### `/respond` -- 템플릿 응답 생성

일반적인 문의 유형에 대해 설정된 템플릿 기반 응답을 생성합니다.

```
/respond [inquiry-type]
```

지원 문의 유형에는 정보주체 요청, discovery hold, 벤더 질문, NDA 요청, 직접 정의한 custom category가 포함됩니다.

## 스킬

| 스킬 | 설명 |
|-------|-------------|
| `contract-review` | Playbook 기반 계약 분석, 표준 이탈 분류, redline 생성을 수행합니다. |
| `nda-triage` | NDA 선별 기준, 분류 규칙, 라우팅 추천을 제공합니다. |
| `compliance` | 개인정보 규제(GDPR, CCPA), DPA 검토, 정보주체 요청을 다룹니다. |
| `canned-responses` | 템플릿 관리, 응답 범주, escalation trigger를 다룹니다. |
| `legal-risk-assessment` | escalation criteria가 포함된 severity-by-likelihood framework로 법무 위험을 평가하고 분류합니다. 계약 위험 평가, 거래 exposure 평가, 이슈 심각도 분류, matter가 senior counsel 또는 outside legal review를 필요로 하는지 판단할 때 사용합니다. |
| `meeting-briefing` | 법무 관련성이 있는 미팅을 위한 구조화된 브리핑을 준비하고 그에 따른 action item을 추적합니다. 계약 협상, 이사회 미팅, 컴플라이언스 검토처럼 법무 맥락, 배경 조사, action tracking이 필요한 미팅 준비에 사용합니다. |

## 예시 워크플로

### 계약 검토

1. 이메일로 벤더 계약을 받습니다.
2. `/review-contract`를 실행하고 문서를 업로드합니다.
3. 맥락을 제공합니다: "우리는 고객이고, 분기 말까지 close해야 하며 data protection과 liability에 집중해야 합니다."
4. GREEN/YELLOW/RED flag가 포함된 조항별 분석을 받습니다.
5. YELLOW 및 RED 항목에 대한 구체적인 redline 문구를 받습니다.
6. 분석 결과를 deal team과 공유합니다.

### NDA 분류

1. Sales team이 신규 prospect의 NDA를 보냅니다.
2. `/triage-nda`를 실행하고 NDA를 붙여 넣거나 업로드합니다.
3. 즉시 분류 결과를 받습니다: GREEN(서명 라우팅), YELLOW(검토할 특정 이슈), RED(전체 법무 검토 필요).
4. GREEN NDA는 바로 승인하고, YELLOW/RED는 표시된 이슈를 처리합니다.

### 일일 브리프

1. `/brief daily`로 아침을 시작합니다.
2. 밤사이 들어온 계약 요청, 컴플라이언스 질문, 다가오는 마감일, 법무 준비가 필요한 캘린더 항목 요약을 받습니다.
3. 긴급도와 마감일 기준으로 하루를 우선순위화합니다.

### 벤더 확인

1. Business team이 기존 벤더와의 새 engagement에 대해 묻습니다.
2. `/vendor-check Acme Corp`를 실행합니다.
3. 기존 계약, 만료일, 핵심 조건을 한눈에 확인합니다.
4. 새 NDA가 필요한지, 기존 조건으로 진행할 수 있는지 바로 판단합니다.

## MCP 통합

> 익숙하지 않은 placeholder가 보이거나 연결된 도구를 확인해야 한다면 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

플러그인은 MCP(Model Context Protocol) 서버를 통해 도구에 연결됩니다.

| 범주 | 예시 | 목적 |
|----------|----------|---------|
| Chat | Slack, Teams | 팀 요청, 알림, 분류 |
| Cloud storage | Box, Egnyte | Playbook, template, precedent |
| Office suite | Microsoft 365 | 이메일, 캘린더, 문서 |
| Project tracker | Atlassian (Jira/Confluence) | Matter tracking, task |

[CONNECTORS.md](CONNECTORS.md)에서 CLM, CRM, e-signature 및 추가 선택지를 포함한 지원 통합 전체 목록을 확인하세요.

`.mcp.json`에서 connection을 설정합니다. Tool을 사용할 수 없으면 플러그인은 가능한 범위로 기능을 축소하고 gap을 표시하며 수동 확인을 제안합니다.

## 커스터마이징

### Playbook 설정

Playbook은 계약 검토 시스템의 핵심입니다. `legal.local.md`에 position을 정의하세요.

- **Standard positions**: 선호하는 계약 조건
- **Acceptable ranges**: escalation 없이 합의할 수 있는 범위
- **Escalation triggers**: senior review 또는 outside counsel이 필요한 조건

### Response template

일반 문의용 template을 정의합니다. Template은 variable substitution을 지원하며, 템플릿 응답을 사용하면 안 되는 상황을 위한 내장 escalation trigger를 포함합니다.

### Risk framework

조직의 위험 선호도와 분류 체계에 맞게 위험 평가 matrix를 커스터마이즈합니다.

## 파일 구조

```
legal/
├── .claude-plugin/plugin.json
├── .mcp.json
├── README.md
├── commands/
│   ├── review-contract.md
│   ├── triage-nda.md
│   ├── vendor-check.md
│   ├── brief.md
│   └── respond.md
└── skills/
    ├── contract-review/SKILL.md
    ├── nda-triage/SKILL.md
    ├── compliance/SKILL.md
    ├── canned-responses/SKILL.md
    ├── legal-risk-assessment/SKILL.md
    └── meeting-briefing/SKILL.md
```
