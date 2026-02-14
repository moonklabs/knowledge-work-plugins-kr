# Enterprise Search

Anthropic의 에이전트 데스크톱 애플리케이션인 [Cowork](https://claude.com/product/cowork)를 위해 설계된 엔터프라이즈 검색 플러그인(plugin)입니다. Claude Code에서도 사용할 수 있습니다. 이메일, 채팅, 문서, 위키 등 회사의 모든 도구를 한곳에서 앱 전환 없이 검색할 수 있습니다.

---

## 작동 방식

하나의 검색 쿼리(search query)로 연결된 모든 도구를 동시에 검색합니다. Claude가 질문을 분해하고 각 데이터 소스(data source)에 맞는 검색을 실행한 후, 출처를 포함한 하나의 일관된 답변으로 결과를 통합합니다.

```
사용자: "API 재설계에 대해 무엇을 결정했나요?"
              ↓ Claude가 모든 소스 검색
~~chat: 결론이 포함된 지난 화요일 #engineering 스레드
~~email: 사양이 포함된 Sarah의 후속 이메일
~~cloud storage: 업데이트된 API 설계 문서 (어제 수정됨)
              ↓ Claude가 결과 통합
"팀은 지난 화요일에 GraphQL 대신 REST를 사용하기로 결정했습니다.
 Sarah님이 목요일에 업데이트된 사양을 보냈고, 설계 문서에도
 최종 접근 방식이 반영되어 있습니다."
```

탭 전환도 필요 없고 어떤 도구에 무엇이 있는지 기억할 필요도 없습니다. 질문하면 답을 얻을 수 있습니다.

---

## 검색 대상

> 익숙하지 않은 플레이스홀더가 보이거나 연결된 도구를 확인하려면 [CONNECTORS.md](CONNECTORS.md)를 참조하세요.

데이터 소스를 자유롭게 조합하여 연결할 수 있습니다. 더 많이 연결할수록 더 완전한 답변을 얻을 수 있습니다.

| 소스 | 검색 내용 |
|------|-----------|
| **~~chat** | 메시지, 스레드, 채널, DM |
| **~~email** | 이메일, 첨부파일, 대화 |
| **~~cloud storage** | 문서, 스프레드시트, 슬라이드, PDF |
| **위키 / 지식 베이스** | 내부 문서, 런북 |
| **프로젝트 관리** | 작업, 이슈, 에픽, 마일스톤 |
| **CRM** | 계정, 연락처, 기회 |
| **티켓팅** | 지원 티켓, 고객 이슈 |

각 소스는 MCP 연결입니다. MCP 설정에서 더 많은 소스를 추가하여 Claude의 검색 범위를 확장할 수 있습니다.

---

## 커맨드

| 커맨드 | 기능 |
|--------|------|
| `/search` | 연결된 모든 소스를 하나의 검색 쿼리로 검색 |
| `/digest` | 모든 소스의 활동에 대한 일간 또는 주간 다이제스트 생성 |

### Search

```
/enterprise-search:search Project Aurora의 상태는 어떤가요?
/enterprise-search:search from:sarah about:budget after:2025-01-01
/enterprise-search:search 이번 주 #product 채널에서 내린 결정 사항
```

`from:`, `in:`, `after:`, `before:`, `type:` 필터를 지원하며, 각 소스의 기본 쿼리 구문에 맞게 지능적으로 적용됩니다.

### Digest

```
/enterprise-search:digest --daily      # 오늘 하루 모든 소스에서 발생한 일
/enterprise-search:digest --weekly     # 프로젝트/주제별로 그룹화된 주간 요약
```

실행 항목, 의사결정 사항, 사용자에 대한 멘션을 강조 표시합니다. 활동을 주제별로 그룹화하여 중요한 내용을 빠르게 확인할 수 있습니다.

---

## 스킬

세 가지 스킬이 검색 경험을 구동합니다:

**검색 전략(Search Strategy)** — 쿼리 분해 및 소스별 변환을 담당합니다. 자연어 질문을 소스별 맞춤 검색으로 분해하고, 모호성을 처리하며, 소스를 사용할 수 없을 때 적절히 대체 전략을 실행합니다.

**소스 관리(Source Management)** — 사용 가능한 MCP 소스를 파악하고, 새로운 소스 연결을 안내하며, 소스 우선순위를 관리하고, 속도 제한을 처리합니다.

**지식 통합(Knowledge Synthesis)** — 여러 소스의 검색 결과를 일관된 답변으로 결합합니다. 소스 간 중복 정보를 제거하고, 출처를 표시하며, 최신성과 신뢰도에 기반한 확신도를 평가하고, 대규모 검색 결과를 요약합니다.

---

## 워크플로우 예시

### 의사결정 찾기

```
사용자: /enterprise-search:search Postgres로 전환하기로 언제 결정했나요?

Claude 검색 중:
  ~~chat → #engineering, #infrastructure 채널에서 "postgres" "전환" "결정" 검색
  ~~email → 제목에 "postgres"가 포함된 스레드 검색
  ~~cloud storage → 데이터베이스 마이그레이션이 언급된 문서 검색

결과: "결정은 3월 3일에 #infrastructure 채널에서 내려졌습니다 (링크).
      Sarah님의 3월 4일 이메일에서 일정이 확정되었으며,
      마이그레이션 계획 문서도 3월 5일에 업데이트되었습니다."
```

### 휴가 후 업무 파악

```
사용자: /enterprise-search:digest --weekly

Claude 스캔 중:
  ~~chat → 내가 속한 채널, DM, 멘션
  ~~email → 받은 편지함 활동
  ~~cloud storage → 공유받거나 수정된 문서

결과: 실행 항목(action items) 및 의사결정 사항이 강조 표시된
      프로젝트별 주간 요약 보고서 생성
```

### 전문가 찾기

```
사용자: /enterprise-search:search 우리 Kubernetes 설정에 대해 누가 잘 아나요?

Claude 검색 중:
  ~~chat → Kubernetes, k8s, 클러스터 관련 메시지
  ~~cloud storage → 인프라 관련 문서 작성자 확인
  위키 → 런북 및 아키텍처 문서 확인

결과: "메시지 이력과 문서 작성 기록을 바탕으로 볼 때,
      Alex와 Priya가 k8s 전문가입니다.
      관심 있는 내용은 이 런북(링크)를 참고하세요."
```

---

## 시작하기

```bash
# 1. 설치
claude plugins add knowledge-work-plugins/enterprise-search

# 2. 모든 소스 통합 검색
/enterprise-search:search [이곳에 질문 입력]

# 3. 요약 리포트 받기
/enterprise-search:digest --daily
```

MCP를 통해 더 많은 소스를 연결할수록 검색 결과가 더 완전해집니다. ~~chat, ~~email, ~~cloud storage부터 시작한 후, 필요에 따라 위키, 프로젝트 관리 도구, CRM을 추가하세요.

---

## 철학

지식 근로자는 매주 여러 도구에 흩어진 정보를 찾는 데 많은 시간을 소비합니다. 답은 어딘가에 존재합니다 — Slack 스레드, 이메일 체인, 문서, 위키 페이지 등 — 하지만 이를 찾으려면 각 도구를 개별적으로 검색하고 결과를 교차 참조하며 올바른 곳을 확인했기를 바라야 합니다.

Enterprise Search는 모든 도구를 하나의 검색 가능한 지식 베이스로 취급합니다. 하나의 쿼리, 모든 소스, 통합된 결과. 회사의 지식이 사일로에 갇혀 있어서는 안 됩니다. 모든 것을 한 번에 검색하세요.
