# Product Management 플러그인

[Cowork](https://claude.com/product/cowork)를 위한 프로덕트 매니지먼트 플러그인입니다. Anthropic의 에이전트형 데스크톱 앱인 Cowork에 맞춰 설계했지만 Claude Code에서도 동작합니다. 기능 스펙 작성, 로드맵 관리, 이해관계자 커뮤니케이션, 사용자 리서치 합성, 경쟁사 분석, 제품 지표 추적까지 PM 워크플로 전반을 다룹니다.

## 설치

```
claude plugins add knowledge-work-plugins/product-management
```

## 제공 기능

이 플러그인은 AI 기반 프로덕트 매니지먼트 파트너로서 다음을 지원합니다.

- **기능 스펙 & PRD**: 문제 정의 또는 기능 아이디어에서 구조화된 PRD를 생성합니다. 사용자 스토리, 요구사항 우선순위, 성공 지표, 스코프 관리가 포함됩니다.
- **로드맵 계획**: 로드맵을 생성, 업데이트, 재우선순위화합니다. Now/Next/Later, 분기 테마, OKR 연계 포맷과 의존성 매핑을 지원합니다.
- **이해관계자 업데이트**: 대상(임원, 엔지니어링, 고객)에 맞춘 상태 업데이트를 생성합니다. 연결된 도구에서 컨텍스트를 가져와 주간 업데이트 부담을 줄입니다.
- **사용자 리서치 합성**: 인터뷰 노트, 설문 데이터, 지원 티켓을 구조화된 인사이트로 변환합니다. 테마를 식별하고 페르소나를 구축하며 근거와 함께 기회 영역을 도출합니다.
- **경쟁 분석**: 경쟁사를 조사하고 기능 비교, 포지셔닝 분석, 전략적 시사점을 담은 브리프를 생성합니다.
- **지표 리뷰**: 제품 지표를 분석하고 추세를 파악하며 목표 대비 비교와 실행 가능한 인사이트를 제공합니다.

## 커맨드

| Command | 설명 |
|---|---|
| `/write-spec` | 문제 정의로부터 기능 스펙 또는 PRD 작성 |
| `/roadmap-update` | 로드맵 생성/업데이트/재우선순위화 |
| `/stakeholder-update` | 이해관계자 업데이트(주간, 월간, 런칭) 생성 |
| `/synthesize-research` | 인터뷰, 설문, 티켓을 기반으로 리서치 합성 |
| `/competitive-brief` | 경쟁 분석 브리프 생성 |
| `/metrics-review` | 제품 지표 리뷰 및 분석 |

## 스킬

| Skill | 범위 |
|---|---|
| `feature-spec` | PRD 구조, 사용자 스토리, 요구사항 분류, 승인 기준 |
| `roadmap-management` | 우선순위 프레임워크(RICE, MoSCoW), 로드맵 포맷, 의존성 매핑 |
| `stakeholder-comms` | 대상별 업데이트 템플릿, 리스크 커뮤니케이션, 의사결정 기록 |
| `user-research-synthesis` | 테마 분석, 어피니티 매핑, 페르소나 개발, 기회 규모화 |
| `competitive-analysis` | 기능 비교 매트릭스, 포지셔닝 분석, win/loss 분석 |
| `metrics-tracking` | 제품 지표 계층, 목표 설정(OKRs), 대시보드 설계, 리뷰 주기 |

## 예시 워크플로

### PRD 작성

```
You: /write-spec
Claude: 어떤 기능이나 문제에 대한 스펙인가요?
You: 엔터프라이즈 고객을 위한 SSO 지원을 추가해야 합니다
Claude: [타겟 사용자, 제약사항, 성공 지표를 질문]
Claude: [문제 정의, 사용자 스토리, 요구사항, 성공 지표, 오픈 이슈가 포함된 PRD 생성]
```

### 이해관계자 업데이트 준비

```
You: /stakeholder-update
Claude: 어떤 유형의 업데이트인가요? (weekly, monthly, launch, ad-hoc)
You: 임원 대상 주간 업데이트입니다
Claude: [프로젝트 관리 도구, 채팅, 문서에서 컨텍스트 수집]
Claude: [진척, 결정 사항, 리스크, 다음 마일스톤이 포함된 요약 생성]
```

### 사용자 리서치 합성

```
You: /synthesize-research
Claude: 어떤 리서치를 합성할까요? 인터뷰 노트를 붙여넣거나 파일을 업로드하거나 연결된 소스에서 가져올 수 있습니다.
You: [인터뷰 녹취 8건 붙여넣기]
Claude: [인터뷰 전반의 테마, 패턴, 인사이트를 식별]
Claude: [핵심 발견, 페르소나, 기회 영역, 근거를 포함한 합성 결과 생성]
```

### 경쟁 분석

```
You: /competitive-brief
Claude: 어떤 경쟁사(또는 기능 영역)를 분석할까요?
You: 우리 온보딩 플로우를 [competitor A]와 [competitor B]와 비교해 주세요
Claude: [경쟁사의 온보딩 접근을 조사]
Claude: [기능 비교, 강점/약점, 전략적 시사점을 포함한 브리프 생성]
```

## 데이터 소스

> 낯선 플레이스홀더가 보이거나 어떤 도구가 연결되어 있는지 확인하려면 [CONNECTORS.md](CONNECTORS.md)를 보세요.

최상의 경험을 위해 프로젝트 관리 및 커뮤니케이션 도구를 연결하세요. 연결하지 않으면 컨텍스트를 수동으로 제공해야 합니다.

**기본 포함 MCP 연결:**
- 채팅(Slack): 팀 컨텍스트 및 이해관계자 스레드
- 프로젝트 트래커(Linear, Asana, monday.com, ClickUp, Atlassian): 로드맵 연동, 티켓 컨텍스트, 상태 추적
- 지식베이스(Notion): 기존 스펙, 리서치, 미팅 노트
- 디자인(Figma): 디자인 컨텍스트 및 핸드오프
- 제품 분석(Amplitude, Pendo): 사용 데이터, 지표, 행동 분석
- 사용자 피드백(Intercom): 지원 티켓, 기능 요청, 사용자 대화
- 미팅 녹취(Fireflies): 미팅 노트 및 논의 컨텍스트

**추가 옵션:**
- 카테고리별 대체 도구는 [CONNECTORS.md](CONNECTORS.md)에서 확인하세요
