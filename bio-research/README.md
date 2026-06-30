# Bio-Research 플러그인

전임상 연구 도구와 데이터베이스(문헌 검색, 유전체 분석, 타깃 우선순위화)에 연결해 초기 단계 생명과학 R&D를 가속합니다. [Cowork](https://claude.com/product/cowork)와 함께 사용하거나 Claude Code에 직접 설치할 수 있습니다.

이 플러그인은 생명과학 연구자를 위해 11개의 MCP 서버 연동과 5개의 분석 스킬을 하나의 패키지로 묶습니다.

## 포함 항목

### MCP 서버(데이터 소스 및 도구)

> 익숙하지 않은 placeholder가 보이거나 연결된 도구를 확인해야 한다면 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

| 제공자 | 수행 내용 | Category/Placeholder |
|----------|-------------|---------------------|
| U.S. National Library of Medicine | 생의학 문헌과 연구 논문 검색 | `~~literature` |
| deepsense.ai | bioRxiv 및 medRxiv preprint 접근 | `~~literature` |
| Consensus | 동료심사 연구의 AI 기반 검색 및 종합 | `~~literature` |
| John Wiley & Sons | 학술 연구 및 출판물 접근 | `~~journal access` |
| Sage Bionetworks | 공동 연구 데이터 관리 | `~~data repository` |
| deepsense.ai | Bioactive drug-like compound 데이터베이스 | `~~chemical database` |
| OpenTargets | Drug target discovery 및 prioritization | `~~drug targets` |
| deepsense.ai | NIH/NLM 임상시험 registry | `~~clinical trials` |
| BioRender | Scientific illustration 생성 | `~~scientific illustration` |
| Owkin | Biology용 AI: histopathology 및 drug discovery | `~~AI research` |
| Benchling\* | 실험실 데이터 관리 플랫폼 | `~~lab platform` |

### 선택형 바이너리 MCP 서버

다음 항목은 별도 바이너리 다운로드가 필요합니다.

- **10X Genomics txg-mcp** (`~~genomics platform`) — 클라우드 분석 데이터 및 작업 흐름 ([GitHub](https://github.com/10XGenomics/txg-mcp/releases))
- **ToolUniverse** (`~~tool database`) — Harvard MIMS의 scientific discovery용 AI tool ([GitHub](https://github.com/mims-harvard/ToolUniverse/releases))

### 스킬(분석 작업 흐름)

#### Single-Cell RNA QC
scverse best practice를 따르는 scRNA-seq data 자동 품질 관리입니다. MAD 기반 필터링과 종합 시각화를 포함하며 `.h5ad`, `.h5` file을 지원합니다.

#### scvi-tools
Single-cell omics용 deep learning toolkit입니다. 통합, 배치 보정, 라벨 전이, 멀티모달 분석을 위해 scVI, scANVI, totalVI, PeakVI, MultiVI, DestVI, veloVI, sysVI model을 다룹니다.

#### Nextflow Pipelines
로컬 또는 공개 GEO/SRA sequencing data에서 nf-core bioinformatics pipeline을 실행합니다.
- **rnaseq** — 유전자 발현 및 차등 발현
- **sarek** — Germline 및 somatic variant calling(WGS/WES)
- **atacseq** — 염색질 접근성 분석

#### Instrument Data to Allotrope
실험실 장비 출력 파일(PDF, CSV, Excel, TXT)을 Allotrope Simple Model(ASM) 형식으로 변환합니다. Cell counter, spectrophotometer, plate reader, qPCR, chromatography system 등 40개 이상의 장비 유형을 지원합니다.

#### Scientific Problem Selection
Fischbach & Walsh framework 기반으로 연구 문제를 선택하기 위한 체계적 framework입니다. 아이디어 발상, 위험평가, 최적화, 의사결정 트리, 역경 계획, 종합을 다루는 9개 skill을 포함합니다.

## 시작하기

```bash
# plugin 설치
/install anthropics/knowledge-work-plugins bio-research

# 사용 가능한 도구 확인
/start
```

## 일반 작업 흐름

**문헌 검토**
~~literature database에서 논문을 검색하고, ~~journal access로 전문에 접근하며, ~~scientific illustration으로 figure를 만듭니다.

**Single-Cell Analysis**
scRNA-seq data에 QC를 실행한 뒤 scvi-tools로 통합, 배치 보정, 세포 유형 주석을 수행합니다.

**Sequencing Pipeline**
GEO/SRA에서 공개 데이터를 다운로드하고 nf-core pipeline(RNA-seq, variant calling, ATAC-seq)을 실행한 뒤 출력을 검증합니다.

**Drug Discovery**
~~chemical database에서 생리활성 화합물을 검색하고, ~~drug target database로 표적 우선순위화를 수행하며 임상시험 데이터를 검토합니다.

**Research Strategy**
Scientific problem selection framework를 사용해 새 아이디어를 pitch하고, 막힌 project를 troubleshoot하거나 전략적 결정을 평가합니다.

## 라이선스

Skill은 Apache 2.0 license를 따릅니다. MCP server는 각 author가 제공하므로 조건은 개별 server documentation을 참고하세요.
