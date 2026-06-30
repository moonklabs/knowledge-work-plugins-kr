# PDF Viewer 플러그인

실시간 인터랙티브 뷰어에서 PDF를 보고, 주석을 달고, 서명합니다. 계약서를 마크업하고,
시각적 피드백으로 양식을 채우고, 승인 스탬프를 찍고,
서명을 배치한 뒤 주석이 적용된 사본을 다운로드합니다.

## 주요 기능

- **PDF 열기** — local file 또는 academic source(arXiv, bioRxiv 등)에서 PDF를 엽니다.
- **Collaborative annotation** — Claude가 highlight, note, stamp를 section별로 제안하고, 사용자가 viewer에서 batch마다 review합니다.
- **Form 작성** — live preview와 함께 field-by-field로 안내하며 작성합니다.
- **Document sign** — page에 signature 또는 initials image를 배치합니다.
- **Approval stamp** — APPROVED, DRAFT, CONFIDENTIAL 또는 custom label을 찍습니다.
- **Download** — viewer toolbar에서 annotated PDF를 export합니다.

## 명령

| 명령 | 수행 내용 |
|---------|-------------|
| `/pdf-viewer:open` | Interactive viewer에서 PDF를 엽니다. |
| `/pdf-viewer:annotate` | Document를 함께 훑으며 markup을 제안, 적용, review합니다. |
| `/pdf-viewer:fill-form` | PDF form field를 interactive하게 채웁니다. |
| `/pdf-viewer:sign` | Page에 signature 또는 initials image를 배치합니다. |

## 일반 PDF 읽기와의 차이

이 플러그인은 문서를 직접 보고 마크업한 뒤 주석이 적용된 사본을 다운로드하려는 **인터랙티브 시각 작업 흐름**용입니다.

PDF에서 **summary 또는 text extraction**만 원한다면 이 plugin을 사용하지 않아도 됩니다. Claude는 PDF file을 native로 읽을 수 있고, pure ingestion에는 그 방식이 더 빠릅니다.

## 작동 방식

이 plugin은 `npx`를 통해 사용자의 machine에서 실행되는 **local MCP server**(`@modelcontextprotocol/server-pdf`)를 사용합니다. API key나 remote service가 필요 없으며, plugin이 load될 때 PDF server가 자동으로 시작됩니다.

## 요구 사항

- Node.js >= 18
- Remote PDF(arXiv 등)를 위한 internet

## 지원 PDF source

- Local file(working directory 안의 file path)
- [arXiv](https://arxiv.org) — `/abs/` URL은 PDF로 auto-convert됩니다
- Direct HTTPS PDF URL(bioRxiv, Zenodo, OSF 등). Landing page가 아니라 PDF link를 사용하세요

## Signature disclaimer

`/pdf-viewer:sign`은 page에 **visual** signature image를 배치합니다. Certified 또는 cryptographic digital signature가 아닙니다. Legally binding e-signature에는 dedicated signing service를 사용하세요.
