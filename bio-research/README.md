# Bio-Research 플러그인

Preclinical research tool과 database(literature search, genomics analysis, target prioritization)에 연결해 early-stage life sciences R&D를 가속합니다. [Cowork](https://claude.com/product/cowork)와 함께 사용하거나 Claude Code에 직접 설치할 수 있습니다.

이 plugin은 life science researcher를 위해 11개의 MCP server integration과 5개의 analysis skill을 하나의 package로 묶습니다.

## 포함 항목

### MCP server(data source 및 tool)

> 익숙하지 않은 placeholder가 보이거나 연결된 tool을 확인해야 한다면 [CONNECTORS.md](CONNECTORS.md)를 참고하세요.

| Provider | 수행 내용 | Category/Placeholder |
|----------|-------------|---------------------|
| U.S. National Library of Medicine | Biomedical literature와 research article 검색 | `~~literature` |
| deepsense.ai | bioRxiv 및 medRxiv preprint 접근 | `~~literature` |
| Consensus | Peer-reviewed research의 AI-powered search 및 synthesis | `~~literature` |
| John Wiley & Sons | Academic research 및 publication 접근 | `~~journal access` |
| Sage Bionetworks | Collaborative research data management | `~~data repository` |
| deepsense.ai | Bioactive drug-like compound database | `~~chemical database` |
| OpenTargets | Drug target discovery 및 prioritization | `~~drug targets` |
| deepsense.ai | NIH/NLM clinical trial registry | `~~clinical trials` |
| BioRender | Scientific illustration 생성 | `~~scientific illustration` |
| Owkin | Biology용 AI — histopathology 및 drug discovery | `~~AI research` |
| Benchling\* | Lab data management platform | `~~lab platform` |

### 선택형 binary MCP server

다음 항목은 별도 binary download가 필요합니다.

- **10X Genomics txg-mcp** (`~~genomics platform`) — Cloud analysis data 및 workflow ([GitHub](https://github.com/10XGenomics/txg-mcp/releases))
- **ToolUniverse** (`~~tool database`) — Harvard MIMS의 scientific discovery용 AI tool ([GitHub](https://github.com/mims-harvard/ToolUniverse/releases))

### 스킬 (Analysis Workflows)

#### Single-Cell RNA QC
scverse best practice를 따르는 scRNA-seq data automated quality control입니다. MAD-based filtering과 comprehensive visualization을 포함해 `.h5ad`, `.h5` file을 지원합니다.

#### scvi-tools
Single-cell omics용 deep learning toolkit입니다. Integration, batch correction, label transfer, multi-modal analysis를 위해 scVI, scANVI, totalVI, PeakVI, MultiVI, DestVI, veloVI, sysVI model을 다룹니다.

#### Nextflow Pipelines
Local 또는 public GEO/SRA sequencing data에서 nf-core bioinformatics pipeline을 실행합니다.
- **rnaseq** — Gene expression 및 differential expression
- **sarek** — Germline 및 somatic variant calling(WGS/WES)
- **atacseq** — Chromatin accessibility analysis

#### Instrument Data to Allotrope
Laboratory instrument output file(PDF, CSV, Excel, TXT)을 Allotrope Simple Model(ASM) format으로 변환합니다. Cell counter, spectrophotometer, plate reader, qPCR, chromatography system 등 40개 이상 instrument type을 지원합니다.

#### Scientific Problem Selection
Fischbach & Walsh framework 기반 research problem selection을 위한 systematic framework입니다. Ideation, risk assessment, optimization, decision tree, adversity planning, synthesis를 다루는 9개 skill을 포함합니다.

## 시작하기

```bash
# plugin 설치
/install anthropics/knowledge-work-plugins bio-research

# 사용 가능한 tool 확인
/start
```

## 일반 workflow

**Literature Review**
~~literature database에서 paper를 검색하고, ~~journal access로 full-text에 접근하며, ~~scientific illustration으로 figure를 만듭니다.

**Single-Cell Analysis**
scRNA-seq data에 QC를 실행한 뒤 scvi-tools로 integration, batch correction, cell type annotation을 수행합니다.

**Sequencing Pipeline**
GEO/SRA에서 public data를 download하고 nf-core pipeline(RNA-seq, variant calling, ATAC-seq)을 실행한 뒤 output을 verify합니다.

**Drug Discovery**
~~chemical database에서 bioactive compound를 검색하고, ~~drug target database로 target prioritization을 수행하며 clinical trial data를 review합니다.

**Research Strategy**
Scientific problem selection framework를 사용해 new idea를 pitch하고, 막힌 project를 troubleshoot하거나 strategic decision을 evaluate합니다.

## 라이선스

Skill은 Apache 2.0으로 licensed됩니다. MCP server는 각 author가 제공하므로 terms는 individual server documentation을 참고하세요.
