# Common Room 플러그인

Common Room 기반 GTM 워크플로입니다. 계정 조사, 연락처 조사, 통화 준비, 개인화된 아웃리치, 잠재 고객 발굴, 주간 브리핑을 지원합니다.

## 개요

이 플러그인은 Claude를 Common Room MCP 서버에 연결하고, 가장 일반적인 영업 담당자 작업 흐름을 다루는 여섯 가지 스킬을 제공합니다. 모든 출력은 실제 Common Room 신호 데이터에 근거합니다. 자사 제품 신호, 파트너 커뮤니티 신호, 외부 의도 신호, RoomieAI와 Spark의 보강 데이터를 활용합니다.

## 요구 사항

- **Common Room MCP**(`mcp.commonroom.io/mcp`)가 연결되고 인증되어 있어야 합니다. 모든 플러그인 기능의 기본 데이터 소스입니다.
- **Calendar 커넥터**(선택) — `call-prep`와 `weekly-prep-brief`에서 자동 미팅 조회를 가능하게 합니다. 연결하지 않으면 두 스킬 모두 미팅 세부 정보를 사용자에게 묻습니다.

## 스킬

스킬은 대화로 트리거됩니다. 원하는 일을 설명하면 Claude가 적절한 스킬을 자동으로 불러옵니다.

| 스킬 | 트리거 문구 |
|-------|----------------|
| `account-research` | Common Room 데이터를 사용해 회사를 조사합니다. 'research [company]', 'tell me about [domain]', 'pull up signals for [account]', 'what's going on with [company]' 또는 계정 수준 질문에서 트리거됩니다. |
| `contact-research` | Common Room 데이터를 사용해 특정 사람을 조사합니다. 'who is [name]', 'look up [email]', 'research [contact]', 'is [name] a warm lead' 또는 연락처 수준 질문에서 트리거됩니다. |
| `call-prep` | Common Room 신호를 사용해 고객 또는 잠재 고객 통화를 준비합니다. 'prep me for my call with [company]', 'prepare for a meeting with [company]', 'what should I know before talking to [company]' 또는 통화 준비 요청에서 트리거됩니다. |
| `compose-outreach` | Common Room 신호를 사용해 개인화된 아웃리치 메시지를 생성합니다. 'draft outreach to [person]', 'write an email to [name]', 'compose a message for [contact]' 또는 아웃리치 초안 작성 요청에서 트리거됩니다. |
| `prospect` | Common Room Prospector를 사용해 타깃 계정 또는 연락처 목록을 만듭니다. 'find companies that match [criteria]', 'build a prospect list', 'find contacts at [type of company]', 'show me companies hiring [role]' 또는 목록 작성 요청에서 트리거됩니다. |
| `weekly-prep-brief` | 다음 7일 동안의 모든 외부 통화에 대한 종합 주간 브리핑을 생성합니다. 'weekly prep brief', 'prepare my week', 'what calls do I have this week', 'Monday prep' 또는 주간 계획 요청에서 트리거됩니다. |

## 명령

명시적 호출이 더 적합한 복잡한 워크플로용 명령 두 가지입니다.

| 명령 | 사용법 |
|---------|-------|
| `/generate-account-plan <company>` | 이해관계자 매핑, 참여 분석, 기회, 위험, 액션 아이템이 포함된 종합 전략 계정 계획 |
| `/weekly-brief [date range]` | 전체 주간 준비 브리핑을 생성합니다. 기본값은 다음 7일입니다 |

## 각 스킬의 출력

**계정 조사** — 전체 개요, 표적 필드 질문, 데이터가 부족할 때의 정직한 응답, MCP 데이터와 LLM 추론 결합 등 네 가지 패턴을 처리합니다. 최신 뉴스용 웹 검색을 포함하며 "My Segments"로 자동 범위를 잡습니다.

**연락처 조사** — 이메일, 이름+회사, 소셜 핸들로 조회합니다. 보강된 신원, CRM 필드, 점수, 웹사이트 방문, 활동 이력, Spark 분석, 대화 시작점을 반환합니다.

**통화 준비** — 회사 스냅샷, 참석자별 프로필, 신호 하이라이트, 맞춤형 대화 포인트, 예상 반론, 추천 통화 결과를 제공합니다. Gong/통화 녹음 활동을 우선하며, 연결된 경우 캘린더를 고려해 동작합니다.

**아웃리치 작성** — Common Room 신호와 웹 검색 hook에 근거한 세 가지 개인화 형식(이메일, 통화 스크립트, LinkedIn 메시지)을 제공합니다. 가능하면 사용자의 회사 포지셔닝에 맞춥니다.

**잠재 고객 발굴** — 신규 회사(ProspectorOrganization)와 기존 계정(Organization)을 구분합니다. 반복 개선과 "find companies like [X]" 같은 유사 회사 검색을 지원합니다. 웹 검색이 신규 결과를 보강합니다.

**주간 준비 브리프** — 다음 7일의 모든 외부 통화를 다루는 전체 브리핑입니다. 회사 스냅샷, 참석자 프로필, 신호, 미팅별 추천 목표를 포함합니다.

## 설정

1. Cowork 설정에서 Common Room MCP 서버가 연결 및 인증된 상태인지 확인합니다.
2. 선택 사항으로 통화 준비 및 주간 브리핑의 자동 미팅 조회를 위해 calendar MCP 서버를 연결합니다.
3. 이 플러그인을 설치합니다. 모든 스킬과 명령을 즉시 사용할 수 있습니다.

## 사용자 컨텍스트

사용자 담당 영역으로 범위를 잡는 모든 스킬은 Common Room에서 `Me` object를 자동으로 가져옵니다. 이를 통해 사용자 프로필, 역할, "My Segments"를 얻고 query가 기본적으로 해당 담당 영역을 사용하게 합니다. 자세한 내용은 `references/me-context.md`를 참고하세요.

회사 컨텍스트가 있으면 스킬이 사용자 제품과 ICP에 맞게 추천을 조정합니다. 자세한 내용은 `references/my-company-context.md`를 참고하세요.

## 커스터마이징

Calendar 커넥터와 도구 참조 동작 방식은 `CONNECTORS.md`를 참고하세요.
