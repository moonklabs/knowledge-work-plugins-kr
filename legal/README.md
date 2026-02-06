# Legal Productivity 플러그인

[Cowork](https://claude.com/product/cowork)를 위한 사내 법무팀용 AI 생산성 플러그인입니다. Anthropic의 에이전트형 데스크톱 앱인 Cowork에 맞춰 설계했지만 Claude Code에서도 동작합니다. 계약 검토, NDA 분류, 컴플라이언스 워크플로, 법무 브리핑, 템플릿 응답을 자동화하며, 조직의 플레이북과 리스크 허용 범위에 맞게 구성할 수 있습니다.

> **면책 고지:** 이 플러그인은 법무 업무를 보조하지만 법률 자문을 제공하지 않습니다. 결론은 반드시 자격을 갖춘 법률 전문가와 검증해야 합니다. AI가 생성한 분석은 법적 결정에 의존하기 전에 변호사가 검토해야 합니다.

## 대상 페르소나

- **Commercial Counsel**: 계약 협상, 벤더 관리, 딜 지원
- **Product Counsel**: 제품 리뷰, 서비스 약관, 개인정보 보호 정책, IP 이슈
- **Privacy / Compliance**: 데이터 보호 규정, DPA 검토, 정보주체 요청, 규제 모니터링
- **Litigation Support**: 디스커버리 홀드, 문서 리뷰 준비, 사건 브리핑

## 설치

```
claude plugins add knowledge-work-plugins/legal
```

## 빠른 시작

### 1. 플러그인 설치

```
claude plugins add knowledge-work-plugins/legal
```

### 2. 플레이북 구성

조직의 표준 입장을 정의하는 로컬 설정 파일을 만듭니다. 여기에는 팀의 협상 플레이북, 리스크 허용 범위, 표준 조항을 기록합니다.

프로젝트의 `.claude/` 디렉터리에 `legal.local.md` 파일을 생성하세요.

```markdown
# 법무 플레이북 구성

## 계약 검토 기준

### 책임 제한
- 표준 입장: 상호 12개월치 지급/지급예정 수수료 한도
- 허용 범위: 6-24개월치 수수료
- 에스컬레이션 트리거: 무제한 책임, 결과적 손해 포함

### 면책(Indemnification)
- 표준 입장: IP 침해와 데이터 침해에 대한 상호 면책
- 허용 범위: 제3자 청구에 한정된 면책
- 에스컬레이션 트리거: 일방 면책 의무, 무제한 면책

### IP 소유권
- 표준 입장: 각 당사자는 기존 IP를 유지, 고객 데이터는 고객 소유
- 에스컬레이션 트리거: 광범위한 IP 양도 조항, 기존 IP에 대한 업무제공(Work-for-Hire) 조항

### 데이터 보호
- 표준 입장: 개인 데이터 처리 시 DPA 필수
- 요구사항: 하위 처리자 통지, 종료 시 데이터 삭제, 72시간 내 침해 통지
- 에스컬레이션 트리거: DPA 미제공, 보호장치 없는 국외 이전

### 기간 및 종료
- 표준 입장: 연 단위 계약 + 편의 해지 30일 통지
- 허용 범위: 최초 기간 이후 편의 해지가 가능한 다년 계약
- 에스컬레이션 트리거: 통지 기간 없는 자동 갱신, 편의 해지 불가

### 준거법
- 선호: [Your jurisdiction]
- 허용: 주요 상사 관할 (NY, DE, CA, England & Wales)
- 에스컬레이션 트리거: 비표준 관할, 불리한 지역에서의 강제 중재

## NDA 기본값
- 상호 의무 필수
- 기간: 기본 2-3년, 영업비밀은 5년
- 표준 제외: 독자 개발, 공개 정보, 제3자로부터 정당하게 수령
- 잔존(Residuals) 조항: 좁은 범위라면 허용

## 응답 템플릿
공통 문의에 대해 템플릿 파일 경로를 지정하거나 인라인 템플릿을 정의하세요.
```

### 3. 도구 연결

플러그인은 MCP로 기존 도구와 연결할 때 가장 효과적입니다. 기본 구성 서버에는 Slack, Box, Egnyte, Atlassian, Microsoft 365가 포함됩니다. 지원되는 카테고리와 옵션의 전체 목록은 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

## 커맨드

### `/review-contract` -- 플레이북 기반 계약 검토

조직의 협상 플레이북에 따라 계약을 검토합니다. 편차를 표시하고, 레드라인을 생성하며, 비즈니스 영향 분석을 제공합니다.

```
/review-contract
```

입력 방식: 파일 업로드, URL, 계약 텍스트 붙여넣기. 컨텍스트(당사/상대방, 마감, 집중 영역)를 질문한 뒤 조항별로 플레이북과 비교합니다.

### `/triage-nda` -- NDA 사전 분류

수신 NDA를 표준 기준으로 빠르게 분류합니다. GREEN(표준 승인), YELLOW(법무 검토), RED(중대한 이슈)로 구분합니다.

```
/triage-nda
```

### `/vendor-check` -- 벤더 계약 상태 확인

연결된 시스템 전반에서 특정 벤더와의 기존 계약 상태를 확인합니다.

```
/vendor-check [vendor name]
```

기존 NDA, MSA, DPA, 만료일, 핵심 조항을 보고합니다.

### `/brief` -- 법무팀 브리핑

법무 업무를 위한 컨텍스트 브리핑을 생성합니다.

```
/brief daily          # 법무 관련 아침 브리핑
/brief topic [query]  # 특정 법무 질문에 대한 리서치 브리핑
/brief incident       # 진행 중인 사안에 대한 긴급 브리핑
```

### `/respond` -- 템플릿 응답 생성

설정된 템플릿에서 공통 문의 유형에 대한 응답을 생성합니다.

```
/respond [inquiry-type]
```

지원 유형: 정보주체 요청, 디스커버리 홀드, 벤더 문의, NDA 요청, 그리고 사용자가 정의한 커스텀 카테고리.

## 스킬

| Skill | 설명 |
|-------|-------------|
| `contract-review` | 플레이북 기반 계약 분석, 편차 분류, 레드라인 생성 |
| `nda-triage` | NDA 스크리닝 기준, 분류 규칙, 라우팅 추천 |
| `compliance` | 개인정보 규정(GDPR, CCPA), DPA 검토, 정보주체 요청 |
| `canned-responses` | 템플릿 관리, 응답 카테고리, 에스컬레이션 트리거 |
| `legal-risk-assessment` | 리스크 심각도 프레임워크, 분류 레벨, 에스컬레이션 기준 |
| `meeting-briefing` | 미팅 준비 방법론, 컨텍스트 수집, 액션 아이템 추적 |

## 예시 워크플로

### 계약 검토

1. 이메일로 벤더 계약을 수신
2. `/review-contract` 실행 후 문서 업로드
3. 컨텍스트 제공: "우리는 고객 측이며, 분기 말까지 클로징 필요, 데이터 보호와 책임 한도에 집중"
4. GREEN/YELLOW/RED 플래그가 포함된 조항별 분석 수신
5. YELLOW/RED 항목에 대한 구체적 레드라인 문구 제공
6. 딜 팀과 분석 결과 공유

### NDA 분류

1. 세일즈 팀이 신규 잠재 고객의 NDA를 전송
2. `/triage-nda` 실행 후 NDA를 붙여넣거나 업로드
3. 즉시 분류 결과 확인: GREEN(서명 라우팅), YELLOW(특정 이슈 검토), RED(전면 법무 검토 필요)
4. GREEN은 즉시 승인, YELLOW/RED는 표시된 이슈를 처리

### 데일리 브리핑

1. 아침에 `/brief daily` 실행
2. 야간 계약 요청, 컴플라이언스 질문, 마감 일정, 법무 준비가 필요한 캘린더 항목 요약 확인
3. 긴급도와 마감 기준으로 하루 우선순위 정리

### 벤더 체크

1. 사업팀이 기존 벤더와의 신규 협업에 대해 문의
2. `/vendor-check Acme Corp` 실행
3. 기존 계약, 만료일, 핵심 조항을 한눈에 확인
4. 신규 NDA가 필요한지, 기존 조건으로 진행 가능한지 즉시 판단

## MCP 통합

> 낯선 플레이스홀더가 보이거나 어떤 도구가 연결되어 있는지 확인하려면 [CONNECTORS.md](CONNECTORS.md)를 보세요.

플러그인은 MCP(Model Context Protocol) 서버를 통해 도구에 연결합니다.

| 카테고리 | 예시 | 목적 |
|----------|----------|---------|
| Chat | Slack, Teams | 팀 요청, 알림, 분류 |
| Cloud storage | Box, Egnyte | 플레이북, 템플릿, 선례 |
| Office suite | Microsoft 365 | 이메일, 캘린더, 문서 |
| Project tracker | Atlassian (Jira/Confluence) | 사건 추적, 업무 |

CLM, CRM, 전자서명 등 추가 옵션을 포함한 전체 지원 통합 목록은 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

`.mcp.json`에서 연결을 구성합니다. 도구를 사용할 수 없으면 플러그인은 그 차이를 알려 주고 수동 확인을 제안합니다.

## 커스터마이징

### 플레이북 구성

플레이북은 계약 검토 시스템의 핵심입니다. `legal.local.md`에 기준을 정의하세요.

- **표준 입장**: 선호하는 계약 조항
- **허용 범위**: 에스컬레이션 없이 합의 가능한 범위
- **에스컬레이션 트리거**: 시니어 리뷰나 외부 자문이 필요한 조항

### 응답 템플릿

공통 문의 유형에 대한 템플릿을 정의합니다. 템플릿은 변수 치환을 지원하며, 템플릿 응답을 쓰면 안 되는 상황에 대한 내장 에스컬레이션 트리거를 포함합니다.

### 리스크 프레임워크

조직의 리스크 허용도와 분류 체계에 맞게 리스크 평가 매트릭스를 커스터마이즈합니다.

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
