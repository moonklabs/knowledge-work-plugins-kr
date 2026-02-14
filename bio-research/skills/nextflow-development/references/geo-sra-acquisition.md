# GEO/SRA 데이터 획득

NCBI GEO/SRA에서 원시 시퀀싱 데이터를 다운로드하고 nf-core 파이프라인용으로 준비합니다.

**사용 시기:** 발표된 데이터셋을 재분석하거나, 연구 결과를 검증하거나, 공개 코호트와 결과를 비교할 때 사용하십시오.

## 목차

- [워크플로우 개요](#워크플로우-개요)
- [1단계: 연구 정보 조회](#1단계-연구-정보-조회)
- [2단계: 시료 그룹 검토](#2단계-시료-그룹-검토)
- [3단계: FASTQ 파일 다운로드](#3단계-fastq-파일-다운로드)
- [4단계: 샘플시트 생성](#4단계-샘플시트-생성)
- [5단계: nf-core 파이프라인 실행](#5단계-nf-core-파이프라인-실행)
- [지원 파이프라인](#지원-파이프라인)
- [지원 생물종](#지원-생물종)
- [전체 예제](#전체-예제)
- [문제 해결](#문제-해결)

---

## 워크플로우 개요

예시: "GSE309891에서 차등 발현 유전자(differential expression genes) 찾기 (약물 처리 vs 대조군)"

```
┌─────────────────────────────────────────────────────────────────┐
│                    GEO/SRA DATA ACQUISITION                     │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
                 ┌────────────────────────┐
                 │   Fetch study info     │
                 │   • Query NCBI/SRA     │
                 │   • Get metadata       │
                 │   • Detect organism    │
                 │   • Identify data type │
                 └────────────────────────┘
                              │
                              ▼
                 ┌────────────────────────┐
                 │   Present summary      │
                 │   • Organism: Human    │
                 │   • Genome: GRCh38     │
                 │   • Type: RNA-Seq      │
                 │   • Pipeline: rnaseq   │
                 │   • Samples: 12        │
                 │     (6 treated,        │
                 │      6 control)        │
                 │   • Size: ~24 GB       │
                 └────────────────────────┘
                              │
                              ▼
                    ┌─────────────────┐
                    │  USER CONFIRMS  │◄──── Decision point
                    │  genome/pipeline│
                    └─────────────────┘
                              │
                              ▼
                 ┌────────────────────────┐
                 │   Select samples       │
                 │   • Group by condition │
                 │   • Show treated/ctrl  │
                 └────────────────────────┘
                              │
                              ▼
                    ┌─────────────────┐
                    │  USER SELECTS   │◄──── Decision point
                    │  sample subset  │
                    └─────────────────┘
                              │
                              ▼
                 ┌────────────────────────┐
                 │   Download FASTQs      │
                 │   • 24 files (R1+R2)   │
                 │   • Parallel transfers │
                 │   • Auto-resume        │
                 └────────────────────────┘
                              │
                              ▼
                 ┌────────────────────────┐
                 │   Generate samplesheet │
                 │   • Map SRR to files   │
                 │   • Pair R1/R2         │
                 │   • Assign conditions  │
                 └────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                    NF-CORE PIPELINE EXECUTION                   │
│              (Continue with Step 1 of main workflow)            │
└─────────────────────────────────────────────────────────────────┘
```

---

## Claude를 위한 지침

GEO/SRA 데이터 획득을 지원할 때:

1. **항상 연구 정보를 먼저 조회하여** 사용자에게 사용 가능한 데이터를 보여주십시오
2. **다운로드 전에 확인을 요청하십시오** — 시료 그룹과 크기를 제시한 후 AskUserQuestion을 사용하여 어떤 하위집합을 다운로드할지 질문하십시오
3. **생물종과 데이터 유형에 따라 적절한 게놈과 파이프라인을 제안하십시오**
4. **데이터 준비가 완료되면 메인 SKILL.md 워크플로우로 돌아가십시오**

확인 질문 예시:
```
Question: "Which sample group would you like to download?"
Options:
  - "RNA-Seq:PAIRED (42 samples, ~87 GB)"
  - "RNA-Seq:SINGLE (7 samples, ~4.5 GB)"
  - "All samples (49 samples, ~92 GB)"
```

---

## 1단계: 연구 정보 조회

다운로드 전에 GEO 연구에 대한 메타데이터를 가져옵니다.

```bash
python scripts/sra_geo_fetch.py info <GEO_ID>
```

**예시:**
```bash
python scripts/sra_geo_fetch.py info GSE110004
```

**출력 내용:**
- 연구 제목 및 요약
- 생물종 (자동 제안 게놈 포함)
- 시료 및 실행 수
- 데이터 유형 (RNA-Seq, ATAC-seq 등)
- 예상 다운로드 크기
- 제안 nf-core 파이프라인

**정보를 JSON으로 저장:**
```bash
python scripts/sra_geo_fetch.py info GSE110004 -o study_info.json
```

---

## 2단계: 시료 그룹 검토

데이터 유형과 레이아웃별로 구성된 시료 그룹을 확인합니다. 여러 데이터 유형이 혼합된 연구에 유용합니다.

```bash
python scripts/sra_geo_fetch.py groups <GEO_ID>
```

**출력 예시:**
```
Sample Group          Count Layout     GSM Range                    Est. Size
--------------------------------------------------------------------------------
RNA-Seq                  42 PAIRED     GSM2879618...(42 samples)      87.4 GB
RNA-Seq                   7 SINGLE     GSM2976181-GSM2976187           4.5 GB
--------------------------------------------------------------------------------
TOTAL                    49                                           91.9 GB

Available groups for --subset option:
  1. "RNA-Seq:PAIRED" - 42 samples (~87.4 GB)
  2. "RNA-Seq:SINGLE" - 7 samples (~4.5 GB)
```

**개별 실행 목록 확인:**
```bash
python scripts/sra_geo_fetch.py list <GEO_ID>

# Filter by data type
python scripts/sra_geo_fetch.py list GSE110004 --filter "RNA-Seq:PAIRED"
```

**결정 지점:** 시료 그룹을 검토하십시오. 연구에 여러 데이터 유형이 있는 경우 어떤 하위집합을 다운로드할지 결정하십시오.

---

## 3단계: FASTQ 파일 다운로드

ENA에서 FASTQ 파일을 다운로드합니다 (SRA보다 빠릅니다).

```bash
python scripts/sra_geo_fetch.py download <GEO_ID> -o <OUTPUT_DIR>
```

**옵션:**
- `-o, --output`: 출력 디렉터리 (필수)
- `-i, --interactive`: 대화형으로 시료 그룹 선택
- `-s, --subset`: 데이터 유형별 필터 (예: "RNA-Seq:PAIRED")
- `-p, --parallel`: 병렬 다운로드 수 (기본값: 4)
- `-t, --timeout`: 다운로드 타임아웃(초) (기본값: 600)

### 대화형 모드 (권장)

여러 데이터 유형이 있는 연구에서 대화형 시료 선택을 위해 `-i` 플래그를 사용하십시오:

```bash
python scripts/sra_geo_fetch.py download GSE110004 -o ./fastq -i
```

**대화형 출력:**
```
============================================================
  SELECT SAMPLE GROUP TO DOWNLOAD
============================================================

  [1] RNA-Seq (paired)
      Samples: 42
      GSM: GSM2879618...(42 samples)
      Size: ~87.4 GB

  [2] RNA-Seq (single)
      Samples: 7
      GSM: GSM2976181-GSM2976187
      Size: ~4.5 GB

  [0] Download ALL (49 samples)
------------------------------------------------------------

Enter selection (0-2):
```

### 직접 하위집합 선택

또는 하위집합을 직접 지정하십시오:

```bash
# Download only RNA-Seq paired-end data
python scripts/sra_geo_fetch.py download GSE110004 -o ./fastq \
    --subset "RNA-Seq:PAIRED" --parallel 6
```

**참고:** 다운로드는 기존 파일을 자동으로 건너뜁니다. 중단된 다운로드는 명령을 다시 실행하여 재개할 수 있습니다.

---

## 4단계: 샘플시트 생성

nf-core 파이프라인과 호환되는 샘플시트를 생성합니다.

```bash
python scripts/sra_geo_fetch.py samplesheet <GEO_ID> \
    --fastq-dir <FASTQ_DIR> \
    -o samplesheet.csv
```

**옵션:**
- `-f, --fastq-dir`: 다운로드된 FASTQ 파일이 있는 디렉터리 (필수)
- `-o, --output`: 출력 샘플시트 경로 (기본값: samplesheet.csv)
- `-p, --pipeline`: 대상 파이프라인 (미지정 시 자동 감지)

**예시:**
```bash
python scripts/sra_geo_fetch.py samplesheet GSE110004 \
    --fastq-dir ./fastq \
    -o samplesheet.csv
```

**출력:** 스크립트는 다음을 수행합니다:
1. 대상 파이프라인이 요구하는 형식으로 샘플시트 생성
2. 제안 게놈 참조 표시
3. 제안 nf-core 명령어 표시

---

## 5단계: nf-core 파이프라인 실행

샘플시트 생성 후 스크립트가 제안 명령어를 제공합니다.

**출력 예시:**
```
Suggested command:
   nextflow run nf-core/rnaseq \
       --input samplesheet.csv \
       --outdir results \
       --genome R64-1-1 \
       -profile docker
```

**결정 지점:** 검토 및 확인하십시오:
1. 제안된 파이프라인이 올바른가?
2. 게놈 참조가 사용자의 생물종에 맞는가?
3. 추가 파이프라인 옵션이 필요한가?

그런 다음 메인 SKILL.md 워크플로우 (1단계: 환경 확인)로 돌아가 파이프라인 실행을 계속하십시오.

---

## 지원 파이프라인

이 스킬은 라이브러리 전략에 따라 적절한 파이프라인을 자동 감지합니다. 별표(★) 표시된 파이프라인은 구성, 샘플시트 생성 및 문서가 완전히 지원됩니다. 그 외에는 제안만 제공되며 nf-core 문서를 참고하여 수동 설정이 필요합니다.

| 라이브러리 전략 | 제안 파이프라인 | 지원 수준 |
|----------------|----------------|----------|
| RNA-Seq        | nf-core/rnaseq | ★ 완전 지원 |
| ATAC-seq       | nf-core/atacseq | ★ 완전 지원 |
| WGS/WXS        | nf-core/sarek  | ★ 완전 지원 |
| ChIP-seq       | nf-core/chipseq | 수동 설정 |
| Bisulfite-Seq  | nf-core/methylseq | 수동 설정 |
| miRNA-Seq      | nf-core/smrnaseq | 수동 설정 |
| Amplicon       | nf-core/ampliseq | 수동 설정 |

---

## 지원 생물종

자동 제안 게놈이 있는 일반 생물종:

| 생물종 | 게놈 | 비고 |
|--------|------|------|
| Homo sapiens | GRCh38 | 인간 참조 게놈 |
| Mus musculus | GRCm39 | 마우스 참조 게놈 |
| Saccharomyces cerevisiae | R64-1-1 | 효모 S288C |
| Drosophila melanogaster | BDGP6 | 초파리 |
| Caenorhabditis elegans | WBcel235 | 예쁜꼬마선충 |
| Danio rerio | GRCz11 | 제브라피시 |
| Arabidopsis thaliana | TAIR10 | 애기장대 |
| Rattus norvegicus | Rnor_6.0 | 랫트 |

전체 목록은 `scripts/config/genomes.yaml`을 참조하십시오.

---

## 전체 예제

GSE110004 (효모 RNA-seq) 재분석:

```bash
# 1. Get study info and sample groups
python scripts/sra_geo_fetch.py info GSE110004

# 2. Download with interactive selection
python scripts/sra_geo_fetch.py download GSE110004 -o ./fastq -i
# Select option [1] for RNA-Seq paired-end samples

# 3. Generate samplesheet
python scripts/sra_geo_fetch.py samplesheet GSE110004 \
    --fastq-dir ./fastq \
    -o samplesheet.csv

# 4. Run nf-core/rnaseq (continue with main SKILL.md workflow)
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --outdir results \
    --genome R64-1-1 \
    -profile docker
```

### 비대화형 다운로드 대안

```bash
# Review sample groups first
python scripts/sra_geo_fetch.py groups GSE110004

# Download specific subset directly
python scripts/sra_geo_fetch.py download GSE110004 \
    --subset "RNA-Seq:PAIRED" \
    -o ./fastq \
    --parallel 4
```

---

## 문제 해결

### ENA 다운로드 실패
ENA 다운로드가 실패하면 SRA에서 직접 데이터를 가져와야 할 수 있습니다:

```bash
# Create SRA tools environment
conda create -n sra_tools -c bioconda sra-tools

# Download with prefetch + fasterq-dump
conda run -n sra_tools prefetch SRR6357070
conda run -n sra_tools fasterq-dump SRR6357070 -O ./fastq
```

### SRA 실행 없음
일부 GEO 데이터셋에는 원시 시퀀싱 리드가 아닌 처리된 데이터만 있습니다. 확인 방법:
```bash
python scripts/sra_geo_fetch.py info <GEO_ID>
```
"Runs: 0"이면 해당 데이터셋에 SRA 원시 데이터가 없을 수 있습니다.

### SuperSeries 지원
여러 SubSeries를 포함하는 GEO SuperSeries는 자동으로 처리됩니다. 도구는 다음을 수행합니다:
1. GEO ID가 SuperSeries인지 감지
2. 연결된 BioProject 접근 번호 확인
3. BioProject에서 모든 SRA 실행 조회

예시: GSE110004는 BioProject PRJNA432544에 연결된 SuperSeries입니다.

### 게놈이 인식되지 않음
생물종이 게놈 매핑에 없으면 게놈을 수동으로 지정하십시오:
```bash
# Check available iGenomes
python scripts/manage_genomes.py list

# Or provide custom reference files to nf-core
nextflow run nf-core/rnaseq --fasta /path/to/genome.fa --gtf /path/to/genes.gtf
```

---

## 요구 사항

- Python 3.8+
- `requests` 라이브러리 (선택 사항이지만 권장)
- `pyyaml` 라이브러리 (선택 사항, 게놈 구성용)
- NCBI 및 ENA에 대한 네트워크 접근

선택 의존성 설치:
```bash
pip install requests pyyaml
```
