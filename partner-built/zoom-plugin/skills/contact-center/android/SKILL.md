---

name: contact-center/android
description: >-
  Android용 Zoom Contact Center SDK입니다. Native Android chat/video/ZVA/scheduled callback integration, campaign
  mode, service lifecycle, rejoin handling에 사용합니다.
user-invocable: false
triggers:
  - "contact center android"
  - "zcc android"
  - "zoomccinterface android"
  - "zoomccchatservice"
  - "zoomccvideoservice"
  - "releasezoomccservice"
  - "android rejoin"
---

# Zoom Contact Center SDK - Android

Official docs:
- https://developers.zoom.us/docs/contact-center/android/
- https://marketplacefront.zoom.us/sdk/contact/android/index.html

## Quick Links

1. [concepts/sdk-lifecycle.md](concepts/sdk-lifecycle.md)
2. [examples/service-patterns.md](examples/service-patterns.md)
3. [references/android-reference-map.md](references/android-reference-map.md)
4. [troubleshooting/common-issues.md](troubleshooting/common-issues.md)

## SDK Surface Summary

- SDK manager: `ZoomCCInterface`
- Channel services:
- `getZoomCCChatService()`
- `getZoomCCVideoService()`
- `getZoomCCZVAService()`
- `getZoomCCScheduledCallbackService()`
- Campaign support via web campaign service and campaign metadata.

## Hard Guardrails

- Initialize SDK in `Application.onCreate`.
- Use `ZoomCCItem` to define channel + identifiers.
- Use `entryId` for chat/video/ZVA.
- Use `apiKey` for scheduled callback and campaign mode.
- Release services on teardown.

## Common Chains

- Contact Center app and engagement context: [../../zoom-apps-sdk/SKILL.md](../../zoom-apps-sdk/SKILL.md)
- Contact Center API automation: [../../rest-api/SKILL.md](../../rest-api/SKILL.md)

## Operations

- [RUNBOOK.md](RUNBOOK.md) - 5-minute preflight and debugging checklist.
