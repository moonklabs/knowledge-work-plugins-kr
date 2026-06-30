---

name: build-zoom-meeting-app
description: >-
  Zoom 미팅 흐름을 만들거나 제품 안에 삽입합니다. Meeting SDK 조인, 웹/모바일 미팅 임베드, 미팅 생명주기 흐름 구현, 또는 Meeting SDK와 Video
  SDK 중 선택이 필요할 때 사용합니다.
---

# /build-zoom-meeting-app

Use this skill for embedded meeting experiences and meeting lifecycle implementation.

## Covers

- Meeting SDK selection and platform routing
- Join/auth implementation planning
- Meeting creation plus join flow design
- Web vs native platform considerations
- Meeting SDK vs Video SDK boundary decisions

## Workflow

1. Confirm whether the user wants a Zoom meeting or a custom video session.
2. Route to Meeting SDK if the user needs actual Zoom meetings.
3. Pull in the relevant platform references.
4. Add REST API only for meeting creation, resource management, or reporting.
5. Add webhooks or RTMS only when the use case explicitly needs them.

## Primary References

- [meeting-sdk](../meeting-sdk/SKILL.md)
- [rest-api](../rest-api/SKILL.md)
- [webhooks](../webhooks/SKILL.md)
- [rtms](../rtms/SKILL.md)
- [video-sdk](../video-sdk/SKILL.md)

## Common Mistakes

- Using Video SDK for normal Zoom meeting embeds
- Mixing resource-management APIs into the core join flow without reason
- Skipping platform-specific SDK constraints until too late
