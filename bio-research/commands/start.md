---
description: bio-research 환경을 설정하고 사용 가능한 도구 탐색
---

# Bio-Research 시작

> 익숙하지 않은 플레이스홀더가 보이거나 연결된 도구를 확인해야 할 경우, [CONNECTORS.md](../CONNECTORS.md)를 참조하십시오.

생물학 연구자가 bio-research 플러그인에 익숙해질 수 있도록 안내합니다. 다음 단계를 순서대로 진행하십시오.

## 1단계: 환영 메시지

다음 환영 메시지를 표시합니다:

```
바이오 리서치 플러그인

생명과학을 위한 AI 기반 연구 보조 도구입니다. 이 플러그인은
문헌 검색, 데이터 분석 파이프라인, 과학 전략을
하나의 환경에 통합합니다.
```

## 2단계: 사용 가능한 MCP 서버 확인

사용 가능한 도구를 나열하여 연결된 MCP 서버를 테스트합니다. 결과를 그룹으로 정리합니다:

**문헌 및 데이터 소스:**
- ~~literature 데이터베이스 — 생의학 문헌 검색
- ~~literature 데이터베이스 — 프리프린트 접근 (생물학 및 의학)
- ~~journal access — 학술 출판물
- ~~data repository — 협업 연구 데이터 (Sage Bionetworks)

**약물 발견 및 임상:**
- ~~chemical database — 생리활성 화합물 데이터베이스
- ~~drug target 데이터베이스 — 약물 표적 발견 플랫폼
- ClinicalTrials.gov — 임상시험 레지스트리
- ~~clinical data platform — 임상시험 사이트 순위 및 플랫폼 도움

**시각화 및 AI:**
- ~~scientific illustration — 과학 그림 및 다이어그램 생성
- ~~AI research 플랫폼 — 생물학 AI (병리조직학, 약물 발견)

연결된 서버와 아직 설정되지 않은 서버를 보고합니다.

## 3단계: 사용 가능한 스킬 확인

이 플러그인에서 사용 가능한 분석 스킬을 나열합니다:

| 스킬 | 기능 |
|------|------|
| **단일세포 RNA QC** | MAD 기반 필터링을 사용한 scRNA-seq 데이터 품질 관리 |
| **scvi-tools** | 단일세포 오믹스를 위한 딥러닝 (scVI, scANVI, totalVI, PeakVI 등) |
| **Nextflow 파이프라인** | nf-core 파이프라인 실행 (RNA-seq, WGS/WES, ATAC-seq) |
| **기기 데이터 변환기** | 실험실 기기 출력을 Allotrope ASM 형식으로 변환 |
| **과학적 문제 선정** | 연구 문제 선택을 위한 체계적 프레임워크 |

## 4단계: 선택 사항 설정 — 바이너리 MCP 서버

별도 설치로 사용 가능한 두 개의 추가 MCP 서버를 안내합니다:

- **~~genomics platform** — 클라우드 분석 데이터 및 워크플로우에 접근
  설치: https://github.com/10XGenomics/txg-mcp/releases 에서 `txg-node.mcpb` 다운로드
- **~~tool database** (Harvard MIMS) — 과학 발견을 위한 AI 도구
  설치: https://github.com/mims-harvard/ToolUniverse/releases 에서 `tooluniverse.mcpb` 다운로드

바이너리 파일 다운로드가 필요하며 선택 사항입니다.

## 5단계: 도움 제안

연구자에게 오늘 무엇을 하고 있는지 물어봅니다. 일반적인 워크플로우를 기반으로 시작점을 제안합니다:

1. **문헌 검토** — "~~literature 데이터베이스에서 [주제]에 대한 최근 논문을 검색해 주세요"
2. **시퀀싱 데이터 분석** — "단일세포 데이터에 QC를 실행해 주세요" 또는 "RNA-seq 파이프라인을 설정해 주세요"
3. **약물 발견** — "~~chemical database에서 [단백질]을 표적으로 하는 화합물을 검색해 주세요" 또는 "[질병]에 대한 약물 표적을 찾아 주세요"
4. **데이터 표준화** — "기기 데이터를 Allotrope 형식으로 변환해 주세요"
5. **연구 전략** — "새로운 프로젝트 아이디어를 평가하는 데 도움을 주세요"

사용자의 응답을 기다린 후 적절한 도구와 스킬로 안내합니다.
