---

name: virtual-agent/android
description: >-
  WebView를 통한 Zoom Virtual Agent Android 통합입니다. Java/Kotlin bridge callback, 네이티브 URL 처리,
  support_handoff relay, 생명주기에 안전한 임베딩에 사용합니다.
user-invocable: false
triggers:
  - "virtual agent android"
  - "android webview zva"
  - "zoomCampaignSdk:ready android"
  - "support_handoff android"
  - "javascriptinterface"
---

# Zoom Virtual Agent - Android

Official docs:
- https://developers.zoom.us/docs/virtual-agent/android/

## Quick Links

1. [concepts/webview-lifecycle.md](concepts/webview-lifecycle.md)
2. [examples/js-bridge-patterns.md](examples/js-bridge-patterns.md)
3. [references/android-reference-map.md](references/android-reference-map.md)
4. [troubleshooting/common-issues.md](troubleshooting/common-issues.md)

## Integration Model

- Host campaign URL in Android WebView.
- Inject runtime context (`window.zoomCampaignSdkConfig`).
- Register JavaScript bridge for `exitHandler`, `commonHandler`, `support_handoff`.
- Apply URL policy via `shouldOverrideUrlLoading` and optional multi-window callbacks.

## Hard Guardrails

- Initialize handlers before expecting JS callbacks.
- Treat legacy `openURL` command handling as compatibility path only.
- Prefer DOM links or `window.open` handling plus explicit native routing.

## Chaining

- Product-level patterns: [../SKILL.md](../SKILL.md)
- Contact Center mobile scope: [../../contact-center/android/SKILL.md](../../contact-center/android/SKILL.md)
