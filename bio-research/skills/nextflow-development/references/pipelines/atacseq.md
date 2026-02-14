# nf-core/atacseq

**버전:** 2.1.2

**공식 문서:** https://nf-co.re/atacseq/2.1.2/
**GitHub:** https://github.com/nf-core/atacseq

> **참고:** 새 버전으로 업데이트할 때 [릴리스 페이지](https://github.com/nf-core/atacseq/releases)에서 주요 변경 사항을 확인하고 아래 명령어의 버전을 업데이트하십시오.

## 목차
- [테스트 명령어](#테스트-명령어)
- [샘플시트 형식](#샘플시트-형식)
- [파라미터](#파라미터)
- [출력 파일](#출력-파일)
- [품질 메트릭](#품질-메트릭)

## 테스트 명령어

```bash
nextflow run nf-core/atacseq -r 2.1.2 -profile test,docker --outdir test_atacseq
```

예상 시간: 약 15분, 피크 및 BigWig 트랙을 생성합니다.

## 샘플시트 형식

```csv
sample,fastq_1,fastq_2,replicate
CONTROL,/path/to/ctrl_rep1_R1.fq.gz,/path/to/ctrl_rep1_R2.fq.gz,1
CONTROL,/path/to/ctrl_rep2_R1.fq.gz,/path/to/ctrl_rep2_R2.fq.gz,2
TREATMENT,/path/to/treat_rep1_R1.fq.gz,/path/to/treat_rep1_R2.fq.gz,1
TREATMENT,/path/to/treat_rep2_R1.fq.gz,/path/to/treat_rep2_R2.fq.gz,2
```

| 열 | 필수 | 설명 |
|----|------|------|
| sample | 예 | 조건/그룹 식별자 |
| fastq_1 | 예 | R1 절대 경로 |
| fastq_2 | 예 | R2 절대 경로 (페어드엔드 필수) |
| replicate | 예 | 반복(replicate) 번호 (정수) |

### 차등 분석을 위한 디자인 파일
```csv
sample,condition
CONTROL,control
TREATMENT,treatment
```

`--deseq2_design design.csv`와 함께 사용하십시오.

## 파라미터

### 최소 실행
```bash
nextflow run nf-core/atacseq -r 2.1.2 -profile docker \
    --input samplesheet.csv --outdir results --genome GRCh38 --read_length 50
```

### 주요 파라미터

| 파라미터 | 기본값 | 설명 |
|---------|--------|------|
| `--genome` | - | `GRCh38`, `GRCh37`, `mm10` |
| `--read_length` | 50 | MACS2 최적화를 위한 리드 길이 |
| `--narrow_peak` | true | 좁은 피크 (브로드의 경우 false) |
| `--mito_name` | chrM | 미토콘드리아 염색체 이름 |
| `--keep_mito` | false | 미토콘드리아 리드 유지 |
| `--min_reps_consensus` | 1 | 합의(consensus) 피크를 위한 최소 반복 수 |

### 차등 접근성 분석
```bash
--deseq2_design design.csv
```

## 출력 파일

```
results/
├── bwa/mergedLibrary/
│   ├── *.mLb.mkD.sorted.bam     # Filtered, deduplicated alignments
│   └── bigwig/
│       └── *.bigWig             # Coverage tracks
├── macs2/narrowPeak/
│   ├── *.narrowPeak             # Peak calls
│   └── consensus/
│       └── consensus_peaks.bed  # Merged peaks across replicates
├── deeptools/
│   ├── plotFingerprint/         # Library complexity
│   └── plotProfile/             # TSS enrichment
├── deseq2/                      # If --deseq2_design provided
└── multiqc/
```

**주요 출력:**
- `*.mLb.mkD.sorted.bam`: 분석 준비 완료 정렬
- `*.narrowPeak`: MACS2 피크 호출 (BED 형식)
- `consensus_peaks.bed`: 반복 간 합의 피크
- `*.bigWig`: 게놈 브라우저 트랙

## 품질 메트릭

| 메트릭 | 양호 | 허용 가능 | 불량 |
|--------|------|----------|------|
| 매핑된 리드 | >80% | 60-80% | <60% |
| 미토콘드리아 | <20% | 20-40% | >40% |
| 중복률 | <30% | 30-50% | >50% |
| FRiP | >30% | 15-30% | <15% |
| TSS 농축 | >6 | 4-6 | <4 |

**단편 크기**: 뉴클레오솜 주기성(nucleosomal periodicity)이 나타나야 합니다 (~50bp 뉴클레오솜 프리, ~200bp 모노-뉴클레오솜).

## 다운스트림 분석

```r
library(ChIPseeker)
library(GenomicRanges)
peaks <- import("consensus_peaks.bed")
peakAnno <- annotatePeak(peaks, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
```

**모티프 분석:**
```bash
findMotifsGenome.pl consensus_peaks.bed hg38 motifs/ -size 200
```

## 문제 해결

**낮은 FRiP**: `plotFingerprint/`에서 라이브러리 복잡도를 확인하십시오. 과도한 전위(over-transposition)를 나타낼 수 있습니다.

**피크 수 부족**: `--macs_qvalue 0.1`로 임계값을 낮추거나 `--narrow_peak false`로 브로드 피크를 사용하십시오.

**높은 중복률**: 저입력 시료에서는 정상이며, 파이프라인이 기본적으로 제거합니다.

## 추가 정보

- **전체 파라미터 목록:** https://nf-co.re/atacseq/2.1.2/parameters/
- **출력 문서:** https://nf-co.re/atacseq/2.1.2/docs/output/
- **사용법 문서:** https://nf-co.re/atacseq/2.1.2/docs/usage/
