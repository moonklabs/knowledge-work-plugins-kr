# Enterprise Search 플러그인

[Cowork](https://claude.com/product/cowork)를 위한 엔터프라이즈 검색 플러그인입니다. Anthropic의 에이전트형 데스크톱 앱인 Cowork에 맞춰 설계했지만 Claude Code에서도 동작합니다. 앱을 옮겨 다니지 않고 이메일, 채팅, 문서, 위키 등 회사의 모든 도구를 한 곳에서 검색합니다.

---

## 동작 방식

하나의 쿼리로 연결된 모든 도구를 동시에 검색합니다. Claude는 질문을 분해해 각 소스에 맞는 검색을 실행하고, 출처를 포함해 하나의 응답으로 통합합니다.

```
You: "API 재설계에 대해 어떤 결정을 내렸나요?"
              ↓ Claude searches
~~chat: #engineering thread from Tuesday with the decision
~~email: Follow-up email from Sarah with the spec
~~cloud storage: Updated API design doc (modified yesterday)
              ↓ Claude synthesizes
"팀은 화요일에 GraphQL 대신 REST로 가기로 결정했습니다.
 Sarah가 목요일에 업데이트된 스펙을 보냈고,
 디자인 문서에는 최종 방향이 반영되어 있습니다."
```

탭을 바꿀 필요도, 어느 도구에 무엇이 있는지 기억할 필요도 없습니다. 질문하면 답을 얻습니다.

---

## 검색 대상

> If you see unfamiliar placeholders or need to check which tools are connected, see [CONNECTORS.md](CONNECTORS.md).

가능한 소스를 자유롭게 조합해 연결하세요. 더 많이 연결할수록 답이 더 완성도 높아집니다.

| Source | 무엇을 찾는가 |
|--------|---------------|
| **~~chat** | 메시지, 스레드, 채널, DM |
| **~~email** | 이메일, 첨부파일, 대화 |
| **~~cloud storage** | 문서, 시트, 슬라이드, PDF |
| **Wiki / Knowledge Base** | 내부 문서, 런북 |
| **Project Management** | 작업, 이슈, 에픽, 마일스톤 |
| **CRM** | 계정, 연락처, 기회 |
| **Ticketing** | 지원 티켓, 고객 이슈 |

각 소스는 MCP 연결입니다. MCP 설정에서 소스를 추가하면 Claude가 검색할 수 있는 범위가 확장됩니다.

---

## 커맨드

| Command | 설명 |
|---------|--------------|
| `/search` | 한 번의 쿼리로 모든 연결 소스 검색 |
| `/digest` | 모든 소스의 활동을 일/주간 다이제스트로 생성 |

### Search

```
/enterprise-search:search what's the status of Project Aurora?
/enterprise-search:search from:sarah about:budget after:2025-01-01
/enterprise-search:search decisions made in #product this week
```

지원 필터: `from:`, `in:`, `after:`, `before:`, `type:` — 각 소스의 쿼리 문법에 맞게 지능적으로 적용됩니다.

### Digest

```
/enterprise-search:digest --daily      # What happened today across all sources
/enterprise-search:digest --weekly     # Weekly rollup grouped by project/topic
```

액션 아이템, 의사결정, 나에 대한 언급을 강조합니다. 활동을 주제별로 묶어 중요한 내용만 빠르게 훑을 수 있습니다.

---

## 스킬

세 가지 스킬이 검색 경험을 구성합니다.

**Search Strategy** — 쿼리 분해 및 소스별 변환. 자연어 질문을 소스별 타겟 검색으로 분해하고, 모호함을 처리하며, 소스가 없을 때는 유연하게 대체합니다.

**Source Management** — 사용 가능한 MCP 소스를 파악하고, 새 소스 연결을 안내하며, 소스 우선순위를 관리하고 레이트 리밋을 처리합니다.

**Knowledge Synthesis** — 여러 소스의 결과를 일관된 답변으로 통합합니다. 소스 간 중복을 제거하고, 출처를 표시하며, 최신성과 권위에 따라 신뢰도를 부여하고, 큰 결과 집합을 요약합니다.

---

## 예시 워크플로

### 결정 사항 찾기

```
You: /enterprise-search:search Postgres로 전환하기로 언제 결정했나요?

Claude 검색:
  ~~chat → #engineering, #infrastructure for "postgres" "switch" "decision"
  ~~email → threads with "postgres" in subject
  ~~cloud storage → docs mentioning database migration

Result: "결정은 3월 3일 #infrastructure에서 내려졌습니다(링크).
         Sarah의 3월 4일 이메일이 일정 확인을 제공했고,
         마이그레이션 계획 문서는 3월 5일에 업데이트되었습니다."
```

### 휴가 이후 캐치업

```
You: /enterprise-search:digest --weekly

Claude 스캔:
  ~~chat → channels you're in, DMs, mentions
  ~~email → inbox activity
  ~~cloud storage → docs shared with you or modified

Result: 프로젝트별로 요약이 묶이며, 액션 아이템은 표시되고
        결정 사항이 강조됩니다.
```

### 전문가 찾기

```
You: /enterprise-search:search 우리 Kubernetes 구성에 대해 잘 아는 사람은?

Claude 검색:
  ~~chat → messages about Kubernetes, k8s, clusters
  ~~cloud storage → docs authored about infrastructure
  Wiki → runbooks and architecture docs

Result: "메시지 히스토리와 문서 작성 기록을 보면
         Alex와 Priya가 k8s 담당자입니다.
         주요 런북은 여기입니다(링크)."
```

---

## 시작하기

```bash
# 1. Install
claude plugins add knowledge-work-plugins/enterprise-search

# 2. 전체 검색
/enterprise-search:search [your question here]

# 3. 다이제스트 받기
/enterprise-search:digest --daily
```

MCP로 연결하는 소스가 많을수록 검색 결과가 더 완성도 높아집니다. ~~chat, ~~email, ~~cloud storage부터 시작하고 필요에 따라 위키, 프로젝트 관리 도구, CRM을 추가하세요.

---

## 철학

지식 노동자는 매주 여러 도구에 흩어진 정보를 찾느라 많은 시간을 씁니다. 답은 어딘가에 존재합니다. Slack 스레드, 이메일 체인, 문서, 위키 페이지 등. 하지만 실제로는 각 도구를 개별적으로 검색하고, 결과를 교차 검증하며, 올바른 곳을 확인했길 바라는 과정이 필요합니다.

Enterprise Search는 모든 도구를 하나의 검색 가능한 지식베이스로 취급합니다. 하나의 쿼리로 모든 소스를 검색하고 결과를 합성합니다. 회사의 지식이 사일로에 갇혀서는 안 됩니다. 한 번에 모두 검색하세요.
