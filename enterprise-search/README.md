# Enterprise Search 플러그인

엔터프라이즈 검색 플러그인입니다. 주로 Anthropic의 agentic 데스크톱 애플리케이션인 [Cowork](https://claude.com/product/cowork)용으로 설계되었지만 Claude Code에서도 동작합니다. 앱을 오가며 전환하지 않고 이메일, 채팅, 문서, 위키 등 회사의 모든 도구를 한곳에서 검색합니다.

---

## 작동 방식

하나의 질의가 연결된 모든 도구를 동시에 검색합니다. Claude는 질문을 분해하고 각 출처에서 표적 검색을 실행한 뒤, 출처 표시가 포함된 하나의 일관된 답변으로 결과를 종합합니다.

```text
사용자: "API redesign에 대해 우리가 무엇을 결정했지?"
              ↓ Claude searches
~~chat: 결정이 담긴 화요일 #engineering thread
~~email: Sarah가 보낸 spec follow-up email
~~cloud storage: Updated API design doc (어제 수정)
              ↓ Claude synthesizes
"Team은 화요일에 GraphQL 대신 REST로 가기로 결정했습니다.
 Sarah가 목요일에 updated spec을 보냈고, design doc은
 final approach를 반영하고 있습니다."
```

Tab을 오갈 필요도, 어떤 도구에 무엇이 있었는지 기억할 필요도 없습니다. 질문하면 답을 얻습니다.

---

## 검색 대상

> 익숙하지 않은 placeholder가 보이거나 연결된 도구를 확인해야 한다면 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

출처를 원하는 조합으로 연결하세요. 더 많이 연결할수록 답변이 완전해집니다.

| 출처 | 찾는 내용 |
|--------|---------------|
| **~~chat** | 메시지, 스레드, 채널, DM |
| **~~email** | 이메일, 첨부파일, 대화 |
| **~~cloud storage** | 문서, 시트, 슬라이드, PDF |
| **Wiki / Knowledge Base** | 내부 문서, runbook |
| **Project Management** | 작업, 이슈, epic, milestone |
| **CRM** | 계정, contact, opportunity |
| **Ticketing** | 지원 티켓, 고객 이슈 |

각 출처는 MCP connection입니다. MCP setting에 출처를 더 추가하면 Claude가 검색할 수 있는 범위가 확장됩니다.

---

## 명령

| 명령 | 수행 내용 |
|---------|--------------|
| `/search` | 하나의 query로 모든 연결된 출처를 검색합니다. |
| `/digest` | 모든 출처 활동의 일간 또는 주간 digest를 생성합니다. |

### Search

```text
/enterprise-search:search Project Aurora 상태가 어떻게 됐지?
/enterprise-search:search from:sarah about:budget after:2025-01-01
/enterprise-search:search 이번 주 #product에서 결정된 사항
```

필터 `from:`, `in:`, `after:`, `before:`, `type:`을 지원하며 각 출처의 native query syntax에 맞게 지능적으로 적용합니다.

### Digest

```text
/enterprise-search:digest --daily      # 오늘 모든 출처에서 있었던 일
/enterprise-search:digest --weekly     # Project/topic별로 묶은 주간 요약
```

액션 아이템, 결정, 멘션을 강조합니다. 활동을 topic별로 묶어 중요한 것만 빠르게 훑을 수 있습니다.

---

## 스킬

세 가지 skill이 검색 경험을 구동합니다.

**Search Strategy** — 질의 분해 및 출처별 변환입니다. 자연어 질문을 출처별 표적 검색으로 나누고 모호성을 처리하며 출처를 사용할 수 없을 때 자연스럽게 대체합니다.

**Source Management** — 사용 가능한 MCP 출처를 파악하고, 새 출처 연결을 안내하며, 출처 우선순위와 rate limit을 관리합니다.

**Knowledge Synthesis** — 여러 출처의 결과를 일관된 답변으로 결합합니다. 출처 간 정보를 중복 제거하고, 출처 표시를 제공하며, 최신성과 권위 기준으로 신뢰도를 산정하고 큰 결과 집합을 요약합니다.

---

## 예시 워크플로

### Decision 찾기

```text
You: /enterprise-search:search Postgres로 전환하기로 언제 결정했지?

Claude searches:
  ~~chat → #engineering, #infrastructure for "postgres" "switch" "decision"
  ~~email → threads with "postgres" in subject
  ~~cloud storage → docs mentioning database migration

Result: "Decision은 3월 3일 #infrastructure에서 내려졌습니다(link).
         Sarah의 3월 4일 email이 timeline을 확인했고,
         migration plan doc은 3월 5일 update되었습니다."
```

### 휴가 후 따라잡기

```text
You: /enterprise-search:digest --weekly

Claude scans:
  ~~chat → 참여 중인 channel, DM, mention
  ~~email → inbox activity
  ~~cloud storage → 공유되었거나 수정된 docs

Result: Project별 grouped summary와 flagged action item,
        highlighted decision을 제공합니다.
```

### Expert 찾기

```text
You: /enterprise-search:search 우리 Kubernetes setup을 누가 잘 알아?

Claude searches:
  ~~chat → Kubernetes, k8s, clusters 관련 message
  ~~cloud storage → infrastructure 관련 docs author
  Wiki → runbook 및 architecture docs

Result: "Message history와 doc authorship 기준으로
         Alex와 Priya가 k8s 담당자입니다.
         Main runbook은 여기 있습니다(link)."
```

---

## 시작하기

```bash
# 1. 설치
claude plugins add knowledge-work-plugins/enterprise-search

# 2. 전체 검색
/enterprise-search:search [your question here]

# 3. Digest 받기
/enterprise-search:digest --daily
```

MCP로 더 많은 source를 연결할수록 search result가 완전해집니다. ~~chat, ~~email, ~~cloud storage부터 시작하고 필요에 따라 wiki, project management tool, CRM을 추가하세요.

---

## 철학

Knowledge worker는 도구 전반에 흩어진 정보를 찾는 데 매주 많은 시간을 씁니다. 답은 Slack thread, email chain, doc, wiki page 어딘가에 있지만, 찾으려면 각 도구를 따로 검색하고 result를 cross-reference하며 올바른 곳을 확인했기를 바라야 합니다.

Enterprise Search는 모든 도구를 하나의 searchable knowledge base처럼 다룹니다. One query, all sources, synthesized results. 회사의 knowledge는 silo에 갇혀 있으면 안 됩니다. 한 번에 모두 검색하세요.
