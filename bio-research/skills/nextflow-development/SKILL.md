---
name: nextflow-development
description: 시퀀싱 데이터에 대해 nf-core 생물정보학 파이프라인(rnaseq, sarek, atacseq)을 실행합니다. RNA-seq, WGS/WES 또는 ATAC-seq 데이터를 분석할 때 사용합니다—로컬 FASTQ 파일 또는 GEO/SRA의 공개 데이터셋 모두 지원합니다. nf-core, Nextflow, FASTQ 분석, 변이 호출(variant calling), 유전자 발현(gene expression), 차등 발현(differential expression), GEO 재분석, GSE/GSM/SRR 접근 번호 또는 샘플시트 생성 관련 요청에서 트리거됩니다.
---

# nf-core 파이프라인 배포

로컬 또는 공개 시퀀싱 데이터에 대해 nf-core 생물정보학 파이프라인을 실행합니다.

**대상 사용자:** 대규모 오믹스(omics) 분석—차등 발현, 변이 호출 또는 크로마틴 접근성 분석—을 수행해야 하는 전문 생물정보학 훈련이 없는 벤치 과학자 및 연구자.

## 워크플로우 체크리스트

```
- [ ] 0단계: 데이터 획득 (GEO/SRA에서 가져오는 경우)
- [ ] 1단계: 환경 확인 (반드시 통과해야 함)
- [ ] 2단계: 파이프라인 선택 (사용자와 확인)
- [ ] 3단계: 테스트 프로필 실행 (반드시 통과해야 함)
- [ ] 4단계: 샘플시트 생성
- [ ] 5단계: 구성 및 실행 (사용자와 게놈 확인)
- [ ] 6단계: 출력 검증
```

---

## 0단계: 데이터 획득 (GEO/SRA인 경우만)

**사용자가 로컬 FASTQ 파일을 가지고 있으면 이 단계를 건너뛰십시오.**

공개 데이터셋의 경우 먼저 GEO/SRA에서 가져옵니다. 전체 워크플로우는 [references/geo-sra-acquisition.md](references/geo-sra-acquisition.md)를 참조하십시오.

**빠른 시작:**

```bash
# 1. Get study info
python scripts/sra_geo_fetch.py info GSE110004

# 2. Download (interactive mode)
python scripts/sra_geo_fetch.py download GSE110004 -o ./fastq -i

# 3. Generate samplesheet
python scripts/sra_geo_fetch.py samplesheet GSE110004 --fastq-dir ./fastq -o samplesheet.csv
```

**결정 지점:** 연구 정보를 가져온 후 사용자에게 확인합니다:
- 다운로드할 시료 하위집합 (여러 데이터 유형이 있는 경우)
- 제안된 게놈 및 파이프라인

그런 다음 1단계로 계속 진행합니다.

---

## 1단계: 환경 확인

**먼저 실행하십시오. 환경이 통과하지 않으면 파이프라인이 실패합니다.**

```bash
python scripts/check_environment.py
```

모든 필수 검사를 통과해야 합니다. 실패하면 수정 지침을 제공합니다:

### Docker 문제

| 문제 | 해결 방법 |
|------|----------|
| 미설치 | https://docs.docker.com/get-docker/ 에서 설치 |
| 권한 거부 | `sudo usermod -aG docker $USER` 후 재로그인 |
| 데몬 미실행 | `sudo systemctl start docker` |

### Nextflow 문제

| 문제 | 해결 방법 |
|------|----------|
| 미설치 | `curl -s https://get.nextflow.io \| bash && mv nextflow ~/bin/` |
| 버전 < 23.04 | `nextflow self-update` |

### Java 문제

| 문제 | 해결 방법 |
|------|----------|
| 미설치 / < 11 | `sudo apt install openjdk-11-jdk` |

**모든 검사가 통과할 때까지 진행하지 마십시오.** HPC/Singularity의 경우 [references/troubleshooting.md](references/troubleshooting.md)를 참조하십시오.

---

## 2단계: 파이프라인 선택

**결정 지점: 진행하기 전에 사용자에게 확인하십시오.**

| 데이터 유형 | 파이프라인 | 버전 | 목표 |
|------------|-----------|------|------|
| RNA-seq | `rnaseq` | 3.22.2 | 유전자 발현 |
| WGS/WES | `sarek` | 3.7.1 | 변이 호출 |
| ATAC-seq | `atacseq` | 2.1.2 | 크로마틴 접근성 |

데이터에서 자동 감지:
```bash
python scripts/detect_data_type.py /path/to/data
```

파이프라인별 상세 정보:
- [references/pipelines/rnaseq.md](references/pipelines/rnaseq.md)
- [references/pipelines/sarek.md](references/pipelines/sarek.md)
- [references/pipelines/atacseq.md](references/pipelines/atacseq.md)

---

## 3단계: 테스트 프로필 실행

**소규모 데이터로 환경을 검증합니다. 실제 데이터 전에 반드시 통과해야 합니다.**

```bash
nextflow run nf-core/<pipeline> -r <version> -profile test,docker --outdir test_output
```

| 파이프라인 | 명령어 |
|-----------|--------|
| rnaseq | `nextflow run nf-core/rnaseq -r 3.22.2 -profile test,docker --outdir test_rnaseq` |
| sarek | `nextflow run nf-core/sarek -r 3.7.1 -profile test,docker --outdir test_sarek` |
| atacseq | `nextflow run nf-core/atacseq -r 2.1.2 -profile test,docker --outdir test_atacseq` |

확인:
```bash
ls test_output/multiqc/multiqc_report.html
grep "Pipeline completed successfully" .nextflow.log
```

테스트가 실패하면 [references/troubleshooting.md](references/troubleshooting.md)를 참조하십시오.

---

## 4단계: 샘플시트 생성

### 자동 생성

```bash
python scripts/generate_samplesheet.py /path/to/data <pipeline> -o samplesheet.csv
```

스크립트 기능:
- FASTQ/BAM/CRAM 파일 탐색
- R1/R2 리드 페어링
- 시료 메타데이터 추론
- 쓰기 전 검증

**sarek의 경우:** 자동 감지되지 않으면 종양/정상 상태를 입력하도록 요청합니다.

### 기존 샘플시트 검증

```bash
python scripts/generate_samplesheet.py --validate samplesheet.csv <pipeline>
```

### 샘플시트 형식

**rnaseq:**
```csv
sample,fastq_1,fastq_2,strandedness
SAMPLE1,/abs/path/R1.fq.gz,/abs/path/R2.fq.gz,auto
```

**sarek:**
```csv
patient,sample,lane,fastq_1,fastq_2,status
patient1,tumor,L001,/abs/path/tumor_R1.fq.gz,/abs/path/tumor_R2.fq.gz,1
patient1,normal,L001,/abs/path/normal_R1.fq.gz,/abs/path/normal_R2.fq.gz,0
```

**atacseq:**
```csv
sample,fastq_1,fastq_2,replicate
CONTROL,/abs/path/ctrl_R1.fq.gz,/abs/path/ctrl_R2.fq.gz,1
```

---

## 5단계: 구성 및 실행

### 5a. 게놈 가용성 확인

```bash
python scripts/manage_genomes.py check <genome>
# If not installed:
python scripts/manage_genomes.py download <genome>
```

일반 게놈: GRCh38 (인간), GRCh37 (레거시), GRCm39 (마우스), R64-1-1 (효모), BDGP6 (초파리)

### 5b. 결정 지점

**결정 지점: 사용자에게 확인하십시오:**

1. **게놈:** 사용할 참조 게놈
2. **파이프라인별 옵션:**
   - **rnaseq:** 정렬기(aligner) (star_salmon 권장, 저메모리의 경우 hisat2)
   - **sarek:** 도구 (생식세포 변이에는 haplotypecaller, 체세포 변이에는 mutect2)
   - **atacseq:** read_length (50, 75, 100 또는 150)

### 5c. 파이프라인 실행

```bash
nextflow run nf-core/<pipeline> \
    -r <version> \
    -profile docker \
    --input samplesheet.csv \
    --outdir results \
    --genome <genome> \
    -resume
```

**주요 플래그:**
- `-r`: 버전 고정
- `-profile docker`: Docker 사용 (HPC의 경우 `singularity`)
- `--genome`: iGenomes 키
- `-resume`: 체크포인트에서 계속

**리소스 제한 (필요한 경우):**
```bash
--max_cpus 8 --max_memory '32.GB' --max_time '24.h'
```

---

## 6단계: 출력 검증

### 완료 확인

```bash
ls results/multiqc/multiqc_report.html
grep "Pipeline completed successfully" .nextflow.log
```

### 파이프라인별 주요 출력

**rnaseq:**
- `results/star_salmon/salmon.merged.gene_counts.tsv` — 유전자 카운트
- `results/star_salmon/salmon.merged.gene_tpm.tsv` — TPM 값

**sarek:**
- `results/variant_calling/*/` — VCF 파일
- `results/preprocessing/recalibrated/` — BAM 파일

**atacseq:**
- `results/macs2/narrowPeak/` — 피크 호출
- `results/bwa/mergedLibrary/bigwig/` — 커버리지 트랙

---

## 빠른 참조

일반적인 종료 코드 및 수정 방법은 [references/troubleshooting.md](references/troubleshooting.md)를 참조하십시오.

### 실패한 실행 재개

```bash
nextflow run nf-core/<pipeline> -resume
```

---

## 참조 자료

- [references/geo-sra-acquisition.md](references/geo-sra-acquisition.md) — 공개 GEO/SRA 데이터 다운로드
- [references/troubleshooting.md](references/troubleshooting.md) — 일반적인 문제 및 해결 방법
- [references/installation.md](references/installation.md) — 환경 설정
- [references/pipelines/rnaseq.md](references/pipelines/rnaseq.md) — RNA-seq 파이프라인 상세 정보
- [references/pipelines/sarek.md](references/pipelines/sarek.md) — 변이 호출 상세 정보
- [references/pipelines/atacseq.md](references/pipelines/atacseq.md) — ATAC-seq 상세 정보

---

## 면책 조항

이 스킬은 자동화된 분석 워크플로우를 위해 nf-core 생물정보학 파이프라인을 Claude Code에 통합하는 방법을 보여주는 프로토타입 예제로 제공됩니다. 현재 구현은 세 가지 파이프라인(rnaseq, sarek, atacseq)을 지원하며, 커뮤니티가 nf-core 파이프라인의 전체 세트로 지원을 확장할 수 있는 기반 역할을 합니다.

교육 및 연구 목적으로 제공되며, 특정 사용 사례에 대한 적절한 검증 없이는 프로덕션 준비가 된 것으로 간주해서는 안 됩니다. 사용자는 컴퓨팅 환경이 파이프라인 요구 사항을 충족하는지 확인하고 분석 결과를 검증할 책임이 있습니다.

Anthropic은 생물정보학 출력의 정확성을 보장하지 않으며, 사용자는 계산 분석 검증을 위한 표준 관행을 따라야 합니다. 이 통합은 nf-core 커뮤니티의 공식 승인을 받았거나 제휴 관계에 있는 것이 아닙니다.

## 저작자 표시

결과를 출판할 때 해당 파이프라인을 인용하십시오. 인용 정보는 각 nf-core 저장소의 CITATIONS.md 파일에서 확인할 수 있습니다 (예: https://github.com/nf-core/rnaseq/blob/3.22.2/CITATIONS.md).

## 라이선스

- **nf-core 파이프라인:** MIT License (https://nf-co.re/about)
- **Nextflow:** Apache License, Version 2.0 (https://www.nextflow.io/about-us.html)
- **NCBI SRA Toolkit:** Public Domain (https://github.com/ncbi/sra-tools/blob/master/LICENSE)
