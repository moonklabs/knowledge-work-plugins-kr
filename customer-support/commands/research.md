---
description: 고객 질문/주제에 대한 다중 소스 리서치 및 출처 포함 답변 생성
argument-hint: "<question or topic>"
---

# Research

> 낯선 플레이스홀더가 보이거나 어떤 도구가 연결되어 있는지 확인하려면 [CONNECTORS.md](../CONNECTORS.md)를 보세요.

고객 질문에 대해 문서/티켓/대화 등 다중 소스를 결합해 답변합니다.

## Usage

```bash
/research <question>
```

## Workflow

1. 질문 범위와 요구 정확도 확인
2. 연결 소스 검색(문서, KB, 과거 티켓, 내부 대화)
3. 상충 정보 정리 및 신뢰도 평가
4. 출처 포함 답변과 불확실성 표시

## Output

- Answer summary
- Source attribution
- Confidence level
- Follow-up questions
