---
description: bio-research 환경 초기화 및 사용 가능한 도구 안내
---

# Bio-Research Start

> 낯선 플레이스홀더가 보이거나 어떤 도구가 연결되어 있는지 확인하려면 [CONNECTORS.md](../CONNECTORS.md)를 보세요.

bio-research 플러그인 사용자를 온보딩합니다.

## Step 1: 환영 및 목표 확인

- 현재 연구 목표(문헌조사, 단일세포 분석, 파이프라인 실행, 타깃 발굴 등)를 확인합니다.

## Step 2: 연결 상태 확인

- `~~literature`, `~~clinical trials`, `~~chemical database`, `~~drug targets`, `~~data repository`, `~~lab platform` 연결 여부를 점검합니다.
- 미연결 소스는 명확히 안내하고 대체 흐름을 제시합니다.

## Step 3: 스킬 선택 안내

목표에 따라 권장 스킬:
- 단일세포 QC: `single-cell-rna-qc`
- 통합/배치보정: `scvi-tools`
- 시퀀싱 워크플로: `nextflow-development`
- 장비 데이터 표준화: `instrument-data-to-allotrope`
- 연구 문제 선정: `scientific-problem-selection`

## Step 4: 첫 실행 예시 제시

- 문헌 검색 쿼리 예시
- 데이터 파일 입력 예시
- 기대 산출물(리포트/시각화/워크플로 결과) 예시

## Step 5: 운영 원칙

- 재현성을 위해 입력/파라미터/버전을 기록
- 가정과 한계를 명시
- 고위험 결론은 실험/전문가 검증 권고
