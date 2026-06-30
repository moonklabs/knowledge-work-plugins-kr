# Create an Asset

**모든 영업팀을 위한 스킬**

몇 분 안에 고객에게 바로 공유할 수 있는 전문적인 영업 자료를 생성합니다. 디자인 스킬은 필요 없습니다.

---

## 주요 기능

이 스킬은 다음 정보를 물어본 뒤 맞춤형 영업 자료를 만듭니다.

1. **잠재고객** — 누구인지, 지금까지 무엇을 논의했는지
2. **대상 독자** — 누가 볼 것인지, 무엇을 중요하게 여기는지
3. **목적** — 무엇을 달성하고 싶은지
4. **형식** — 어떤 방식으로 보여주고 싶은지

그 다음 잠재고객을 조사하고, 문구를 쓰고, 디자인하고, 고객에게 공유할 수 있는 완성도 높은 자료를 빌드합니다.

---

## 지원 형식

| 형식 | 적합한 상황 | 산출물 |
|--------|----------|--------|
| **Interactive Landing Page** | 임원 미팅, 가치 제안 presentation | Demo와 calculator가 포함된 multi-tab page |
| **Deck-Style** | 공식 presentation, 큰 audience | Navigation이 있는 linear slide |
| **One-Pager** | Leave-behind, 빠른 요약 | Single-scroll executive summary |
| **Workflow / Architecture Demo** | 기술 deep-dive, POC proposal | Animated flow가 있는 interactive diagram |

---

## Quick Start

### Option 1: 간단한 prompt

```text
Create an asset for Acme Corp
```

### Option 2: Context 포함

```text
Create an asset for Acme Corp. I met with their VP Engineering
last week - they're struggling with slow release cycles and
want to improve developer productivity. This is for a follow-up
with their technical team.
```

### Option 3: Workflow demo

```text
I want to mock up a workflow showing how a customer would use
our product to automate their invoice processing. The flow is:
invoices come in via email → our AI extracts data → validates
against their ERP → flags exceptions for human review.
```

---

## 생성되는 자료

### Interactive Landing Page

- Tab navigation
- 회사 지표와 조사 내용
- Use case demo
- ROI calculator
- 잠재고객의 brand color를 반영한 전문적인 dark theme

### Deck-Style

- 양사 logo가 들어간 title slide
- Agenda
- Content slide(슬라이드당 메시지 하나)
- Summary 및 next steps
- Speaker note 포함

### One-Pager

- Hero statement
- 핵심 value point 3개
- Proof point
- 명확한 CTA

### Workflow Demo

- 시각적 component node
- Animated data flow
- 단계별 walkthrough
- Play/pause/step control
- 각 단계의 sample data

---

## 진행 방식

```text
1. 사용자가 context를 제공합니다(prospect, audience, purpose)
          ↓
2. Skill이 잠재고객 회사를 조사합니다
          ↓
3. Skill이 3-4개의 확인 질문을 합니다
          ↓
4. 사용자가 방향을 확인합니다
          ↓
5. Skill이 asset을 빌드합니다
          ↓
6. 필요에 따라 반복 수정합니다
```

---

## 자료 공유

Output은 self-contained HTML file입니다. 다음 방식으로 공유할 수 있습니다.

- **Static hosting**: Netlify, Vercel, GitHub Pages 또는 임의 web host에 업로드합니다.
- **Password protect**: 대부분의 host는 간단한 password protection을 제공합니다.
- **Direct share**: HTML file을 이메일로 전송하면 offline에서도 동작합니다.
- **Embed**: 다른 page나 portal에 iframe으로 삽입합니다.

---

## 좋은 결과를 위한 팁

### 풍부한 context 제공

이전 대화, pain point, stakeholder concern을 많이 공유할수록 자료가 더 정확히 맞춤화됩니다.

### Transcript 업로드

Call recording, meeting note, email thread가 있으면 업로드하세요. 스킬이 핵심 quote와 priority를 추출합니다.

### Audience를 구체적으로 지정

"Technical team"도 좋지만, "우리 security model을 평가하는 IT architect"처럼 구체적으로 말하면 더 좋습니다.

### 자유롭게 반복 수정

초안이 딱 맞지 않으면 바꿀 내용을 말하세요. Color, section, messaging, flow 모두 조정할 수 있습니다.

---

## Examples

| 시나리오 | 형식 | 핵심 기능 |
|----------|--------|--------------|
| Discovery 이후 임원 미팅 | Interactive page | ROI calculator, 상대가 말한 priority, case study |
| Technical architecture review | Workflow demo | System diagram, data flow, integration point |
| Board presentation leave-behind | One-pager | Strategic alignment, key metric, 명확한 CTA |
| Large stakeholder meeting | Deck-style | Linear narrative, slide당 한 point, appendix |

---

## FAQ

**Q: 어떤 product/company에도 동작하나요?**
A: 네. 스킬은 이메일 domain에서 판매 중인 대상을 감지하고 그에 맞춰 조사합니다.

**Q: 잠재고객의 brand color는 어떻게 알까요?**
A: 잠재고객 website 또는 brand guideline에서 추출합니다. 이후 수정할 수 있습니다.

**Q: 우리 회사 branding을 대신 쓸 수 있나요?**
A: 네. 첫 빌드 후 brand color를 우리 회사 기준으로 바꿔 달라고 요청하면 됩니다.

**Q: Research가 틀리면 어떻게 하나요?**
A: 틀린 부분을 표시하고 correction을 제공하세요. 스킬이 정확한 정보로 다시 생성합니다.

**Q: PDF로 export할 수 있나요?**
A: 네. Print-optimized version을 요청한 뒤 browser의 print-to-PDF를 사용하세요.

---

## Support

질문이나 feedback이 있다면 알려주세요. 이 스킬은 public sales skills collection의 일부입니다.

---

*슬라이드 디자인보다 영업에 집중하고 싶은 salespeople을 위해 만들었습니다.*
