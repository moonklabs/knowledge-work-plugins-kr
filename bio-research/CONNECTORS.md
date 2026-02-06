# 커넥터

## 도구 레퍼런스 동작 방식

플러그인 파일은 `~~category`를 사용자가 연결한 해당 카테고리의 도구를 가리키는 플레이스홀더로 사용합니다. 예를 들어 `~~literature`는 PubMed, bioRxiv 또는 MCP 서버가 있는 다른 문헌 소스를 의미할 수 있습니다.

플러그인은 **도구에 종속되지 않습니다**. 특정 제품명이 아니라 카테고리(문헌, 임상시험, 화합물 데이터베이스 등)로 워크플로를 설명합니다. `.mcp.json`에는 특정 MCP 서버가 미리 구성되어 있지만, 해당 카테고리의 어떤 MCP 서버라도 사용할 수 있습니다.

## 이 플러그인의 커넥터

| 카테고리 | 플레이스홀더 | 포함 서버 | 기타 옵션 |
|----------|-------------|-----------------|---------------|
| 문헌 | `~~literature` | PubMed, bioRxiv | Google Scholar, Semantic Scholar |
| 과학 일러스트레이션 | `~~scientific illustration` | BioRender | — |
| 임상시험 | `~~clinical trials` | ClinicalTrials.gov | EU Clinical Trials Register |
| 화합물 데이터베이스 | `~~chemical database` | ChEMBL | PubChem, DrugBank |
| 약물 타깃 | `~~drug targets` | Open Targets | UniProt, STRING |
| 데이터 저장소 | `~~data repository` | Synapse | Zenodo, Dryad, Figshare |
| 저널 접근 | `~~journal access` | Wiley Scholar Gateway | Elsevier, Springer Nature |
| AI 연구 | `~~AI research` | Owkin | — |
| 실험실 플랫폼 | `~~lab platform` | Benchling\* | — |

\* 플레이스홀더 — MCP URL이 아직 설정되지 않았습니다
