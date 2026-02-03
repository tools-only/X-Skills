---
title: Model Request Hangs Before Stream Opens
link: stream-hang-timeout
type: debug
created_at: 2026-01-21T13:00:00Z
updated_at: 2026-01-21T13:00:00Z
tags:
  - stream
  - timeout
  - model-request
  - debug
---

# Model Request Hangs Before Stream Opens

## Summary
Requests hang in the model request phase and eventually hit the global 120s timeout. History cleanup is no longer the blocker; the hang occurs before streaming opens and persists even after forcing a non-streaming fallback.

## Evidence (Latest Debug)
- History is clean; no tool calls or dangling tool returns.
- Stream init occurs, but stream never opens.
- Stream watchdog fires at 30s and falls back to non-streaming.
- After fallback, request still times out at 120s with GlobalRequestTimeoutError.

Key log lines:
- `Stream init: node=ModelRequestNode request_id=f3551347 iteration=2 ctx_messages=0 ctx_messages_type=None`
- `Stream request parts: count=1 type=ModelRequest`
- `Stream request part[0]: kind=user-prompt content=gm gm (5 chars)`
- `Stream watchdog timeout after 30.0s; falling back to non-streaming`
- `GlobalRequestTimeoutError: Request exceeded global timeout of 120.0s`

## What We Tried
1. Fixed dangling tool-call cleanup to scan the full history and remove unmatched tool calls.
2. Added detailed debug logging for history, request parts, response parts, and tool returns.
3. Added a stream watchdog timeout to abort streaming after 30s and fall back to non-streaming.

## Conclusion
The hang is inside the model API call (or provider streaming start), not history corruption. The request payload is minimal and valid, but the provider never opens a stream or returns a response.

## Next Step (Not Implemented)
Add a hard request-level watchdog/fallback strategy (e.g., cancel and retry with non-streaming or alternate provider) or verify provider health/configuration.
