---
name: month-heads-up
description: 매월 25일에 실행되어 다음 30일 현금흐름 전망을 보여주고 월말 전에 주의가 필요한 항목을 표시합니다. 선택적 30일 또는 60일 기간을 받습니다.
allowed-tools: Read, WebFetch, Bash
---

Run the month-end heads-up. Pull forward-looking cash data and give the owner a clear "here's what the next 30 days look like" picture with specific things to watch.

Parse arguments:
- `--horizon` (default: `30`) — forecast window in days (`30` or `60`)

## Step 1 — Current cash position

Using the `cash-flow-snapshot` skill workflow:

1. Pull QuickBooks current cash and receivables balance.
2. Pull PayPal settled balance and pending payouts.
3. Combine for total available + incoming cash.

## Step 2 — Upcoming obligations

1. Pull recurring expenses from QuickBooks (payroll, subscriptions, rent/lease) due in the next 30 days.
2. Pull any outstanding invoices past due or due within 14 days.
3. Flag any payment that would push the balance below a comfortable buffer (default: <$2,000 or owner's QB average monthly expense × 0.5).

## Step 3 — Cash-flow forecast

1. Project 30-day net cash: current balance + expected inflows − known obligations.
2. Identify the single tightest week (lowest projected balance).
3. Flag if any week projects negative.

## Step 4 — Two things to watch

Surface no more than two specific, actionable watches:
- Which invoice(s) to chase now
- Which expense(s) to defer or negotiate

Format as:

```
Month-End Heads Up — {current date}
Horizon: next {X} days

Cash today: ${amount}
Projected end-of-period: ${amount}
Tightest week: {date range} — projected ${amount}

TWO THINGS TO WATCH
1. {item} — {why it matters} — suggested action: {action}
2. {item} — {why it matters} — suggested action: {action}
```

## Connector failures

If QuickBooks is unreachable, stop — the cash forecast requires QB as the source of truth. If PayPal is missing, run the forecast from QB-only data and note "PayPal not connected — PayPal receivables excluded from forecast." Same for Stripe/Square if missing.

## Approval gates

- **Never initiate payments or send emails automatically.** Surface the data and actions for the owner to take.
- **Never project revenue that hasn't been confirmed in QB or PayPal.** Use conservative estimates only.

## Output

Present the formatted brief and offer to draft chase emails for any flagged overdue invoices.
