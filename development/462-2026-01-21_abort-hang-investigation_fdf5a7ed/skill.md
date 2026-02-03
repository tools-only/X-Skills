---
title: ESC Abort Causes Subsequent Requests to Hang
link: abort-hang-investigation
type: debug
created_at: 2026-01-21T12:00:00Z
updated_at: 2026-01-21T12:00:00Z
tags:
  - abort
  - hang
  - timeout
  - tool-calls
  - message-history
---

# ESC Abort Causes Subsequent Requests to Hang

## Symptom

User hits ESC to abort mid-request. All subsequent requests in that conversation hang forever (timeout after 120s). Issue persists even after reloading the session.

## Investigation Timeline

### Initial Hypothesis: Broken HTTP Client (WRONG)

We thought cancelling an asyncio task might leave the httpx client connection in a bad state.

**What we did:**
1. Added `invalidate_agent_cache(model, state_manager)` function
2. Called it on abort (`UserAbortError`, `CancelledError`)
3. Called it on timeout (`TimeoutError`)

**Why it didn't work:**
On session reload, it's a fresh process - caches are empty anyway. The agent was being recreated fresh, but requests still hung.

### Second Hypothesis: Debug Logging Reveals Truth

Added detailed logging of message history. Key finding:

```
[DEBUG]   msg[-3]: kind=response, parts=[tool-call, tool-call, tool-call]
[DEBUG]   msg[-2]: kind=request, parts=[retry-prompt, retry-prompt, retry-prompt, user-prompt, user-prompt]
[DEBUG]   msg[-1]: kind=request, parts=[user-prompt:5chars]
```

## Root Cause Found

**The message history has UNANSWERED TOOL CALLS in the middle (not at the end).**

Sequence of events:
1. Model makes 3 tool calls (`msg[-3]`: response with tool-call parts)
2. User hits ESC during tool execution
3. Tool calls never completed - no `tool-return` messages added
4. User sends a new message
5. New user message (`msg[-1]`) becomes the last message
6. Session is saved with corrupted history

**Why `_remove_dangling_tool_calls` doesn't catch this:**

```python
def _remove_dangling_tool_calls(messages, ...):
    while messages:
        last_message = messages[-1]
        has_tool_calls = _message_has_tool_calls(last_message)
        if not has_tool_calls:
            break  # <-- STOPS HERE because last msg is user-prompt
        # ... remove logic
```

The function only removes TRAILING messages with unanswered tool calls. If a user message comes AFTER the dangling tool calls, the cleanup stops.

## The Bug

OpenAI/Anthropic APIs require that every tool call has a corresponding tool result. When tool calls exist without results (even in the middle of history), the API:
- May hang waiting for tool results
- May reject the request
- May behave unpredictably

## Proposed Fix

Need to scan the ENTIRE message history for tool calls without corresponding tool returns, not just trailing messages.

**Algorithm:**
1. Collect all `tool_call_id` values from `response` messages with `tool-call` parts
2. Collect all `tool_call_id` values from `request` messages with `tool-return` parts
3. Any tool call without a matching return = dangling
4. Remove the dangling tool call message OR inject a synthetic "aborted" tool return

## Files Modified During Investigation

- `src/tunacode/core/agents/main.py` - Added debug logging, cache invalidation calls
- `src/tunacode/core/agents/agent_components/agent_config.py` - Added `invalidate_agent_cache()`
- `src/tunacode/core/agents/agent_components/__init__.py` - Export new function
- `tests/unit/core/test_agent_cache_abort.py` - Test for cache invalidation (passes but doesn't fix issue)

## What We Learned

1. The bug is in **message history corruption**, not HTTP client state
2. `_remove_dangling_tool_calls` only handles trailing dangling calls
3. User sending a new message after abort "hides" the dangling calls from cleanup
4. The `retry-prompt` parts in the logs suggest some retry logic ran but didn't properly clean up

## Next Steps

1. Fix `_remove_dangling_tool_calls` to scan entire history
2. Or: Add a "full history validation" pass before sending to API
3. Test by reproducing: make tool calls, ESC abort, send new message, verify cleanup
