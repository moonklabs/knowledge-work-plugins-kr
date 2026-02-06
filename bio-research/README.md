# Bio-Research 플러그인

문헌 검색, 유전체 분석, 타깃 우선순위화 등 전임상 연구 도구와 데이터베이스를 연결해 초기 생명과학 R&D를 가속합니다. [Cowork](https://claude.com/product/cowork)에서 사용하거나 Claude Code에 직접 설치할 수 있습니다.

이 플러그인은 생명과학 연구자를 위해 10개의 MCP 서버 통합과 5개의 분석 스킬을 하나로 묶었습니다.

## 구성 요소

### MCP 서버(데이터 소스 및 도구)

> 낯선 플레이스홀더가 보이거나 어떤 도구가 연결되어 있는지 확인하려면 [CONNECTORS.md](CONNECTORS.md)를 보세요.

| 제공자 | 기능 | 카테고리/플레이스홀더 |
|----------|-------------|---------------------|
| U.S. National Library of Medicine | 생의학 문헌 및 연구 논문 검색 | `~~literature` |
| deepsense.ai | bioRxiv, medRxiv 프리프린트 접근 | `~~literature` |
| John Wiley & Sons | 학술 연구 및 출판물 접근 | `~~journal access` |
| Sage Bionetworks | 협업 연구 데이터 관리 | `~~data repository` |
| deepsense.ai | 생체활성 약물 유사 화합물 데이터베이스 | `~~chemical database` |
| OpenTargets | 약물 타깃 탐색 및 우선순위화 | `~~drug targets` |
| deepsense.ai | NIH/NLM 임상시험 레지스트리 | `~~clinical trials` |
| BioRender | 과학 일러스트레이션 제작 | `~~scientific illustration` |
| Owkin | 생물학용 AI — 조직병리 및 신약 발굴 | `~~AI research` |
| Benchling\* | 실험실 데이터 관리 플랫폼 | `~~lab platform` |

### 선택적 바이너리 MCP 서버

다음은 별도의 바이너리 다운로드가 필요합니다.

- **10X Genomics txg-mcp** (`~~genomics platform`) — 클라우드 분석 데이터 및 워크플로 ([GitHub](https://github.com/10XGenomics/txg-mcp/releases))
- **ToolUniverse** (`~~tool database`) — Harvard MIMS의 과학 발견용 AI 도구 ([GitHub](https://github.com/mims-harvard/ToolUniverse/releases))

### 스킬(분석 워크플로)

#### Single-Cell RNA QC
scverse 모범 사례를 따른 scRNA-seq 데이터 자동 품질 관리. MAD 기반 필터링과 종합 시각화를 포함하며 `.h5ad`, `.h5` 파일을 지원합니다.

#### scvi-tools
단일세포 오믹스용 딥러닝 툴킷. 통합, 배치 보정, 라벨 전이, 멀티모달 분석을 위한 scVI, scANVI, totalVI, PeakVI, MultiVI, DestVI, veloVI, sysVI를 포함합니다.

#### Nextflow Pipelines
로컬 또는 공개 GEO/SRA 시퀀싱 데이터에 대해 nf-core 파이프라인을 실행합니다.
- **rnaseq** — 유전자 발현 및 차등 발현 분석
- **sarek** — 생식세포/체성 변이 호출(WGS/WES)
- **atacseq** — 염색질 접근성 분석

#### Instrument Data to Allotrope
실험 장비 출력 파일(PDF, CSV, Excel, TXT)을 Allotrope Simple Model(ASM) 형식으로 변환합니다. 세포 계수기, 분광광도계, 플레이트 리더, qPCR, 크로마토그래피 시스템 등 40종 이상 장비를 지원합니다.

#### Scientific Problem Selection
Fischbach & Walsh 프레임워크에 기반한 연구 문제 선정 프레임워크입니다. 아이데이션, 리스크 평가, 최적화, 의사결정 트리, 역경 계획, 합성 등 9개 스킬을 포함합니다.

## 시작하기

```bash
# 플러그인 설치
/install anthropics/knowledge-work-plugins bio-research

# 사용 가능한 도구 확인
/start
```

## 공통 워크플로

**문헌 리뷰**
~~literature에서 논문을 검색하고, ~~journal access로 원문에 접근하며, ~~scientific illustration으로 도판을 제작합니다.

**단일세포 분석**
scRNA-seq 데이터에 QC를 수행한 뒤 scvi-tools로 통합, 배치 보정, 세포 타입 주석화를 진행합니다.

**시퀀싱 파이프라인**
GEO/SRA에서 공개 데이터를 다운로드하고 nf-core 파이프라인(RNA-seq, 변이 호출, ATAC-seq)을 실행해 결과를 검증합니다.

**신약 발굴**
~~chemical database에서 생체활성 화합물을 검색하고, ~~drug targets로 타깃을 우선순위화한 뒤 임상시험 데이터를 검토합니다.

**연구 전략**
새 아이디어 제안, 정체된 프로젝트의 문제 해결, 전략적 의사결정 평가를 과학적 문제 선정 프레임워크로 수행합니다.

## 라이선스

스킬은 Apache 2.0 라이선스입니다. MCP 서버는 각 작성자가 제공하므로 약관은 개별 서버 문서를 참고하세요.
