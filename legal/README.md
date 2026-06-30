# Legal Productivity 플러그인

In-house legal team을 위한 AI-powered productivity plugin입니다. 주로 Anthropic의 agentic desktop application인 [Cowork](https://claude.com/product/cowork)용으로 설계되었지만 Claude Code에서도 동작합니다. Contract review, NDA triage, compliance workflow, legal briefing, templated response를 자동화하며, 조직의 playbook과 risk tolerance에 맞게 모두 설정할 수 있습니다.

> **면책 고지:** 이 플러그인은 legal workflow를 보조하지만 법률 자문을 제공하지 않습니다. 결론은 항상 자격을 갖춘 법률 전문가와 검증하세요. AI가 생성한 분석은 법적 의사결정에 사용하기 전에 licensed attorney의 검토를 받아야 합니다. 이 플러그인의 기본 playbook 예시는 미국 법적 입장과 관할(Delaware, New York, California)을 반영합니다. 다른 법체계(EU, UK, Netherlands, Australia 등)에서 운영한다면, 플러그인 분석에 의존하기 전에 `.claude/legal.local.md`의 playbook을 해당 관할의 법적 요구사항, 표준 계약 조건, compliance obligation에 맞게 커스터마이즈해야 합니다.

## 대상 사용자

- **Commercial Counsel** -- Contract negotiation, vendor management, deal support를 담당하는 법무 담당자
- **Product Counsel** -- Product review, terms of service, privacy policy, IP matter를 다루는 법무 담당자
- **Privacy / Compliance** -- Data protection regulation, DPA review, data subject request, regulatory monitoring 담당자
- **Litigation Support** -- Discovery hold, document review prep, case briefing 지원 담당자

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

조직의 standard position을 정의하는 local settings file을 만드세요. 팀의 negotiation playbook, risk tolerance, standard terms를 여기에 기록합니다.

Claude가 찾을 수 있는 위치에 `legal.local.md` file을 만드세요.

- **Cowork**: Folder picker로 Cowork와 공유한 아무 folder에 저장합니다. Plugin이 자동으로 찾습니다.
- **Claude Code**: Project의 `.claude/` directory에 저장합니다.

```markdown
# Legal Playbook Configuration

## Contract Review Positions

### Limitation of Liability
- Standard position: Mutual cap at 12 months of fees paid/payable
- Acceptable range: 6-24 months of fees
- Escalation trigger: Uncapped liability, consequential damages inclusion

### Indemnification
- Standard position: Mutual indemnification for IP infringement and data breach
- Acceptable: Indemnification limited to third-party claims only
- Escalation trigger: Unilateral indemnification obligations, uncapped indemnification

### IP Ownership
- Standard position: Each party retains pre-existing IP; customer owns customer data
- Escalation trigger: Broad IP assignment clauses, work-for-hire provisions for pre-existing IP

### Data Protection
- Standard position: Require DPA for any personal data processing
- Requirements: Sub-processor notification, data deletion on termination, breach notification within 72 hours
- Escalation trigger: No DPA offered, cross-border transfer without safeguards

### Term and Termination
- Standard position: Annual term with 30-day termination for convenience
- Acceptable: Multi-year with termination for convenience after initial term
- Escalation trigger: Auto-renewal without notice period, no termination for convenience

### Governing Law
- Preferred: [Your jurisdiction]
- Acceptable: Major commercial jurisdictions (NY, DE, CA, England & Wales)
- Escalation trigger: Non-standard jurisdictions, mandatory arbitration in unfavorable venue

## NDA Defaults
- Mutual obligations required
- Term: 2-3 years standard, 5 years for trade secrets
- Standard carveouts: independently developed, publicly available, rightfully received from third party
- Residuals clause: acceptable if narrowly scoped

## Response Templates
Configure paths to your template files or define inline templates for common inquiries.
```

### 3. 도구 연결

이 plugin은 MCP를 통해 기존 tool에 연결할 때 가장 잘 동작합니다. Pre-configured server에는 Slack, Box, Egnyte, Atlassian, Microsoft 365가 포함됩니다. 지원 category와 option 전체 목록은 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

## 명령

### `/review-contract` -- Playbook 기반 계약 검토

조직의 negotiation playbook 기준으로 contract를 검토합니다. Deviation을 표시하고 redline을 생성하며 business impact analysis를 제공합니다.

```
/review-contract
```

입력: file upload, URL, pasted contract text. Context(your side, deadline, focus areas)를 물은 뒤 설정된 playbook 기준으로 clause-by-clause review를 수행합니다.

### `/triage-nda` -- NDA 사전 screening

Incoming NDA를 standard criteria 기준으로 빠르게 triage합니다. GREEN(standard approval), YELLOW(counsel review), RED(significant issues)로 분류합니다.

```
/triage-nda
```

### `/vendor-check` -- Vendor agreement 상태

Connected system 전반에서 특정 vendor와의 existing agreement 상태를 확인합니다.

```
/vendor-check [vendor name]
```

Existing NDA, MSA, DPA, expiration date, key term을 보고합니다.

### `/brief` -- Legal team briefing

Legal work를 위한 contextual briefing을 생성합니다.

```
/brief daily          # Morning brief of legal-relevant items
/brief topic [query]  # Research brief on a specific legal question
/brief incident       # Rapid brief on a developing situation
```

### `/respond` -- Template response 생성

Common inquiry type에 대해 configured template 기반 response를 생성합니다.

```
/respond [inquiry-type]
```

지원 inquiry type에는 data subject request, discovery hold, vendor question, NDA request, 직접 정의한 custom category가 포함됩니다.

## 스킬

| 스킬 | 설명 |
|-------|-------------|
| `contract-review` | Playbook 기반 contract analysis, deviation classification, redline generation |
| `nda-triage` | NDA screening criteria, classification rules, routing recommendations |
| `compliance` | Privacy regulation(GDPR, CCPA), DPA review, data subject request |
| `canned-responses` | Template management, response category, escalation trigger |
| `legal-risk-assessment` | Escalation criteria가 포함된 severity-by-likelihood framework로 legal risk를 assess하고 classify합니다. Contract risk evaluation, deal exposure assessment, issue severity classification, matter가 senior counsel 또는 outside legal review를 필요로 하는지 판단할 때 사용합니다. |
| `meeting-briefing` | 법무 관련성이 있는 meeting을 위한 structured briefing을 준비하고 그에 따른 action item을 추적합니다. Contract negotiation, board meeting, compliance review처럼 legal context, background research, action tracking이 필요한 meeting 준비에 사용합니다. |

## 예시 워크플로

### Contract Review

1. Email로 vendor contract를 받습니다.
2. `/review-contract`를 실행하고 document를 업로드합니다.
3. Context를 제공합니다: "우리는 customer이고, quarter end까지 close해야 하며 data protection과 liability에 집중해야 합니다."
4. GREEN/YELLOW/RED flag가 포함된 clause-by-clause analysis를 받습니다.
5. YELLOW 및 RED item에 대한 specific redline language를 받습니다.
6. Analysis를 deal team과 공유합니다.

### NDA Triage

1. Sales team이 new prospect의 NDA를 보냅니다.
2. `/triage-nda`를 실행하고 NDA를 paste 또는 upload합니다.
3. 즉시 classification을 받습니다: GREEN(route for signature), YELLOW(specific issues to review), RED(needs full counsel review).
4. GREEN NDA는 바로 approve하고, YELLOW/RED는 flagged issue를 처리합니다.

### Daily Brief

1. `/brief daily`로 아침을 시작합니다.
2. Overnight contract request, compliance question, upcoming deadline, legal prep이 필요한 calendar item summary를 받습니다.
3. Urgency와 deadline 기준으로 하루를 우선순위화합니다.

### Vendor Check

1. Business team이 existing vendor와의 new engagement에 대해 묻습니다.
2. `/vendor-check Acme Corp`를 실행합니다.
3. Existing agreement, expiration date, key term을 한눈에 확인합니다.
4. New NDA가 필요한지, existing terms로 진행할 수 있는지 바로 판단합니다.

## MCP 통합

> 익숙하지 않은 placeholder가 보이거나 연결된 tool을 확인해야 한다면 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

Plugin은 MCP(Model Context Protocol) server를 통해 tool에 연결됩니다.

| Category | 예시 | 목적 |
|----------|----------|---------|
| Chat | Slack, Teams | Team requests, notifications, triage |
| Cloud storage | Box, Egnyte | Playbooks, templates, precedents |
| Office suite | Microsoft 365 | Email, calendar, documents |
| Project tracker | Atlassian (Jira/Confluence) | Matter tracking, tasks |

[CONNECTORS.md](CONNECTORS.md)에서 CLM, CRM, e-signature 및 additional option을 포함한 supported integration 전체 목록을 확인하세요.

`.mcp.json`에서 connection을 설정합니다. Tool을 사용할 수 없으면 plugin은 graceful degrade하며 gap을 표시하고 manual check를 제안합니다.

## 커스터마이징

### Playbook 설정

Playbook은 contract review system의 핵심입니다. `legal.local.md`에 position을 정의하세요.

- **Standard positions**: Your preferred contract terms
- **Acceptable ranges**: What you can agree to without escalation
- **Escalation triggers**: Terms that require senior review or outside counsel

### Response template

Common inquiry용 template을 정의합니다. Template은 variable substitution을 지원하며, templated response를 사용하면 안 되는 상황을 위한 built-in escalation trigger를 포함합니다.

### Risk framework

조직의 risk appetite와 classification scheme에 맞게 risk assessment matrix를 customize합니다.

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
