# nf-core/sarek

**버전:** 3.7.1

**공식 문서:** https://nf-co.re/sarek/3.7.1/
**GitHub:** https://github.com/nf-core/sarek

> **참고:** 새 버전으로 업데이트할 때 [릴리스 페이지](https://github.com/nf-core/sarek/releases)에서 주요 변경 사항을 확인하고 아래 명령어의 버전을 업데이트하십시오.

## 목차
- [테스트 명령어](#테스트-명령어)
- [샘플시트 형식](#샘플시트-형식)
- [변이 호출 모드](#변이-호출-모드)
- [파라미터](#파라미터)
- [출력 파일](#출력-파일)

## 테스트 명령어

```bash
nextflow run nf-core/sarek -r 3.7.1 -profile test,docker --outdir test_sarek
```

예상 시간: 약 20분, 정렬된 BAM 및 변이 호출 결과를 생성합니다.

## 샘플시트 형식

### FASTQ 기반
```csv
patient,sample,lane,fastq_1,fastq_2
patient1,tumor,L001,/path/to/tumor_L001_R1.fq.gz,/path/to/tumor_L001_R2.fq.gz
patient1,tumor,L002,/path/to/tumor_L002_R1.fq.gz,/path/to/tumor_L002_R2.fq.gz
patient1,normal,L001,/path/to/normal_R1.fq.gz,/path/to/normal_R2.fq.gz
```

### BAM/CRAM 기반
```csv
patient,sample,bam,bai
patient1,tumor,/path/to/tumor.bam,/path/to/tumor.bam.bai
patient1,normal,/path/to/normal.bam,/path/to/normal.bam.bai
```

### 종양/정상 상태 포함
```csv
patient,sample,lane,fastq_1,fastq_2,status
patient1,tumor,L001,tumor_R1.fq.gz,tumor_R2.fq.gz,1
patient1,normal,L001,normal_R1.fq.gz,normal_R2.fq.gz,0
```

`status`: 0 = 정상, 1 = 종양

## 변이 호출 모드

### 생식세포 변이(Germline) (단일 시료)
```bash
nextflow run nf-core/sarek -r 3.7.1 -profile docker \
    --input samplesheet.csv --outdir results --genome GRCh38 \
    --tools haplotypecaller,snpeff
```

### 체세포 변이(Somatic) (종양-정상 쌍)
```bash
nextflow run nf-core/sarek -r 3.7.1 -profile docker \
    --input samplesheet.csv --outdir results --genome GRCh38 \
    --tools mutect2,strelka,snpeff
```

### WES (엑솜)
```bash
nextflow run nf-core/sarek -r 3.7.1 -profile docker \
    --input samplesheet.csv --outdir results --genome GRCh38 \
    --wes --intervals /path/to/targets.bed \
    --tools haplotypecaller,snpeff
```

### 공동 생식세포 변이 호출(Joint germline) (코호트)
```bash
--tools haplotypecaller --joint_germline
```

## 파라미터

### 사용 가능한 도구

**생식세포 변이 호출기(Germline callers):**
- `haplotypecaller`: GATK HaplotypeCaller
- `freebayes`: FreeBayes
- `deepvariant`: DeepVariant (GPU 선택 사항)
- `strelka`: Strelka2 생식세포 변이

**체세포 변이 호출기(Somatic callers):**
- `mutect2`: GATK Mutect2
- `strelka`: Strelka2 체세포 변이
- `manta`: 구조적 변이

**CNV 호출기:**
- `ascat`: 복제수 분석
- `controlfreec`: CNV 감지
- `tiddit`: SV 호출

**주석(Annotation):**
- `snpeff`: 기능 주석
- `vep`: Variant Effect Predictor

### 주요 파라미터

| 파라미터 | 기본값 | 설명 |
|---------|--------|------|
| `--tools` | - | 쉼표로 구분된 도구 목록 |
| `--genome` | - | `GRCh38`, `GRCh37` |
| `--wes` | false | 엑솜 모드 (`--intervals` 필요) |
| `--intervals` | - | 표적 영역 BED 파일 |
| `--joint_germline` | false | 코호트 공동 호출 |
| `--skip_bqsr` | false | 염기 품질 재보정 건너뛰기 |

## 출력 파일

```
results/
├── preprocessing/
│   └── recalibrated/           # Analysis-ready BAMs
│       └── *.recal.bam
├── variant_calling/
│   ├── haplotypecaller/        # Germline VCFs
│   ├── mutect2/                # Somatic VCFs (filtered)
│   └── strelka/
├── annotation/
│   └── snpeff/                 # Annotated VCFs
└── multiqc/
```

## 문제 해결

**BQSR 실패**: 게놈에 사용 가능한 알려진 부위(known sites)를 확인하십시오. 비표준 참조의 경우 `--skip_bqsr`로 건너뛰십시오.

**Mutect2에서 변이 없음**: 샘플시트에서 종양/정상 페어링을 확인하십시오 (`status` 열 확인).

**메모리 부족**: WGS의 경우 `--max_memory '128.GB'`를 사용하십시오.

**DeepVariant GPU**: NVIDIA Docker 런타임이 구성되어 있는지 확인하십시오.

## 추가 정보

- **전체 파라미터 목록:** https://nf-co.re/sarek/3.7.1/parameters/
- **출력 문서:** https://nf-co.re/sarek/3.7.1/docs/output/
- **사용법 문서:** https://nf-co.re/sarek/3.7.1/docs/usage/
