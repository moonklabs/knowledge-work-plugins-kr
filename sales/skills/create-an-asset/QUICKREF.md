# 애셋 생성 — 빠른 참조

## 호출
```
/create-an-asset
/create-an-asset [CompanyName]
"Create an asset for [Company]"
```

---

## 입력 요약

| 입력 | 제공할 정보 |
|------|------------|
| **(a) 잠재 고객** | 기업, 담당자, 거래 단계, 페인 포인트, 녹취록 |
| **(b) 대상** | 임원 / 기술팀 / 운영팀 / 혼합 + 관심 사항 |
| **(c) 목적** | 소개 / 후속 조치 / 심층 분석 / 정렬 / POC / 클로즈 |
| **(d) 형식** | 랜딩 페이지 / 덱 / 원페이저 / 워크플로우 데모 |

---

## 형식 선택

| 필요한 경우... | 선택... |
|---------------|---------|
| 인상적인 멀티탭 경험 | **인터랙티브 랜딩 페이지** |
| 미팅에서 발표할 자료 | **덱 스타일** |
| 남겨줄 간단한 요약 | **원페이저** |
| 시스템 연결 방식의 시각화 | **워크플로우 데모** |

---

## 샘플 프롬프트

**기본:**
```
Create an asset for Acme Corp
```

**맥락 포함:**
```
Create an asset for Acme Corp. They're a manufacturing company
struggling with supply chain visibility. Met with their COO
last week. Need something for the exec team.
```

**워크플로우 데모:**
```
Mock up a workflow for Centric Brands showing how they'd use
our product to monitor contract compliance. Components: our AI,
their Snowflake warehouse, and scanned PDF contracts.
```

---

## 생성 후

| 하고 싶은 것 | 이렇게 말하세요 |
|-------------|----------------|
| 색상 변경 | "Use our brand colors instead" |
| 섹션 추가 | "Add a section on security" |
| 축약 | "Make it more concise" |
| 수정 | "The CEO's name is wrong, it's Jane Smith" |
| PDF 변환 | "Give me a print-friendly version" |

---

## 출력

- 독립형 HTML 파일
- 오프라인 작동
- 어디서든 호스팅 가능 (Netlify, Vercel, GitHub Pages 등)
- 호스팅 제공업체를 통한 비밀번호 보호

---

*맥락 제공 -> 질문 응답 -> 애셋 수령 -> 반복 수정. 이것이 전부입니다.*
