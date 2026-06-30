---
name: signature-request
description: 전자서명를 위해 문서를 준비하고 라우팅합니다. 서명 전 체크리스트를 실행하고 서명 순서를 설정한 뒤 체결으로 보냅니다. 계약가 최종 확정되어 signature 준비가 되었거나, 발송 전 법인명, 첨부자료, 서명란을 검증하거나, 순차/병렬el signer가 있는 envelope을 설정할 때 사용합니다.
argument-hint: "<보낼 document 또는 contract>"
---

# /signature-request -- E-Signature Routing

> If you see unfamiliar placeholders or need to check which tools are connected, see [CONNECTORS.md](../../CONNECTORS.md).

Prepare a document for electronic signature — verify completeness, set signing order, and route for execution.

**Important**: This command assists with legal workflows but does not provide legal advice. Verify documents are in final form before sending for signature.

## Usage

```
/signature-request $ARGUMENTS
```

Prepare for signature: @$1

## Workflow

### Step 1: Accept the Document

Accept the document in any format:
- **File upload**: PDF, DOCX
- **URL**: Link to a document in ~~cloud storage or ~~CLM
- **Reference**: "The Acme Corp MSA we finalized yesterday"

### Step 2: Pre-Signature Checklist

Before routing for signature, verify:

```markdown
## Pre-Signature Checklist

- [ ] Document is in final, agreed form (no open redlines)
- [ ] All exhibits and schedules are attached
- [ ] Correct legal entity names on signature blocks
- [ ] Dates are correct or left blank for execution date
- [ ] Signature blocks match the authorized signers
- [ ] Any required internal approvals have been obtained
- [ ] Document has been reviewed by appropriate counsel
```

### Step 3: Configure Signing

Gather signing details:
- **Signers**: Who needs to sign? (names, emails, titles)
- **Signing order**: Sequential or parallel?
- **Internal approval**: Does anyone need to approve before the counterparty signs?
- **CC recipients**: Who should receive a copy of the executed document?

### Step 4: Route for Signature

**If ~~e-signature is connected:**
- Create the signature envelope/request
- Set signing fields and order
- Add any required initials or date fields
- Send for signature

**If not connected:**
- Generate a signing instruction document
- Provide the document formatted for wet signature or manual e-sign
- List all signers with contact information

## Output

```markdown
## Signature Request: [Document Title]

### Document Details
- **Type**: [MSA / NDA / SOW / Amendment / etc.]
- **Parties**: [Party A] and [Party B]
- **Pages**: [X]

### Pre-Signature Check: [PASS / ISSUES FOUND]
[List any issues that need attention before sending]

### Signing Configuration
| Order | Signer | Email | Role |
|-------|--------|-------|------|
| 1 | [Name] | [email] | [Party A Authorized Signatory] |
| 2 | [Name] | [email] | [Party B Authorized Signatory] |

### CC Recipients
- [Name] — [email]

### Status
[Sent for signature / Ready to send / Issues to resolve first]

### Next Steps
- [What to expect after sending]
- [Expected turnaround time]
- [Follow-up if not signed within X days]
```

## Tips

1. **Check entity names carefully** — The most common signing error is incorrect legal entity names.
2. **Verify authority** — Make sure each signer is authorized to bind their organization.
3. **Keep a copy** — Executed copies should be filed in ~~cloud storage or ~~CLM immediately after execution.
