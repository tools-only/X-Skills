# Research – Streaming Abort Context Loss

**Date:** 2026-01-27 22:46:56
**Owner:** claude
**Phase:** Research
**git_commit:** 8247e85e85f4cead37426260bc819b932472a90a

## Goal

Map out why streaming responses are lost when the user aborts the agent during streaming. Identify where the streaming content lives, how it flows through the system, and why it's not persisted on abort.

## Summary

When a user aborts (Ctrl+C or Escape) during agent streaming, the partial response content is **completely lost**. The streaming content lives only in a UI instance variable (`app.current_stream_text`) which is cleared before cancellation completes. There is no mechanism to capture partial responses and persist them to the session state.

---

## Data Flow

### Normal Flow (Completion)
```
┌─────────────────────────────────────────────────────────────────────────────┐
│ 1. LLM streams tokens → streaming_callback (ui/app.py:348)                │
│    └─→ chunks accumulated in app.current_stream_text (UI variable only)   │
│                                                                              │
│ 2. Stream completes → agent_run.all_messages() returns full history        │
│    └─→ _persist_run_messages() merges into conversation.messages          │
│                                                                              │
│ 3. Finally block (ui/app.py:241-256):                                      │
│    └─→ Content rendered to panel via render_agent_response()              │
│    └─→ Panel written to chat_container                                     │
│    └─→ save_session() persists conversation.messages to disk               │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Cancellation Flow (Abort)
```
┌─────────────────────────────────────────────────────────────────────────────┐
│ 1. User presses Escape → action_cancel_stream (ui/app.py:399)             │
│    ├─→ _streaming_cancelled = True                                         │
│    ├─→ current_stream_text = ""  ← PARTIAL TEXT CLEARED HERE               │
│    └─→ _current_request_task.cancel()                                      │
│                                                                              │
│ 2. asyncio.CancelledError raised → streaming.py:378                        │
│    └─→ Logs "Stream cancelled"                                             │
│    └─→ Re-raises (no partial save)                                         │
│                                                                              │
│ 3. UserAbortError caught in main.py:438-465                                │
│    ├─→ remove_dangling_tool_calls()                                        │
│    ├─→ remove_empty_responses()                                            │
│    ├─→ remove_consecutive_requests()                                       │
│    ├─→ Agent cache invalidated                                             │
│    └─→ Re-raises UserAbortError                                            │
│                                                                              │
│ 4. UI layer catches exception (ui/main.py:75)                              │
│    └─→ Returns silently                                                    │
│                                                                              │
│ RESULT: Partial streaming content is LOST - never persisted anywhere       │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Findings

### Relevant files & why they matter:

| File | Lines | Purpose |
|------|-------|---------|
| `src/tunacode/ui/app.py` | 83, 348-364, 399-407, 241-256 | UI layer: `current_stream_text` accumulation, cancellation, rendering |
| `src/tunacode/core/agents/main.py` | 206-214, 434, 438-465 | Core: message persistence, abort exception handler |
| `src/tunacode/core/agents/agent_components/streaming.py` | 122-384, 350 | Streaming: token delta handling, `_debug_raw_stream_accum` |
| `src/tunacode/core/state.py` | 24-61, 246-275 | State: `SessionState` fields, `save_session()` serialization |
| `src/tunacode/core/agents/resume/sanitize.py` | 253-318, 326-355, 358-404 | Resume: cleanup functions for dangling state |
| `src/tunacode/ui/widgets/chat.py` | 80-150 | Chat container: stream tracking (`start_stream`, `update_stream`, `cancel_stream`) |

### Key Architecture Points:

1. **Layer Separation Issue**: Streaming content lives in UI layer (`ui/app.py:83`), but message persistence happens in core layer (`core/agents/main.py:206-214`). Core layer has no access to partial text during abort.

2. **Success-Only Persistence**: `_persist_run_messages()` (main.py:206) is ONLY called in the success path (line 434). The abort path explicitly does NOT call it.

3. **Debug Accumulator Exists**: `_debug_raw_stream_accum` (state.py:60) captures raw stream content for debugging but is:
   - Only populated in debug mode
   - Not persisted to disk in `save_session()`
   - Never used for user content preservation

4. **Pydantic-AI Message Timing**: `agent_run.all_messages()` only returns complete message history AFTER the run completes. During streaming, the ModelResponse is being constructed internally by pydantic-ai and is not yet available.

---

## Key Patterns / Solutions Found

### Pattern 1: Abort-Cleanup-Reraise
**Location:** `core/agents/main.py:438-465`

The current abort handler prioritizes message integrity over preserving work:
```python
except (UserAbortError, asyncio.CancelledError):
    # DON'T persist agent_run messages - they contain dangling tool calls
    # Just clean up what's already in conversation messages
    remove_dangling_tool_calls(conversation.messages, runtime.tool_registry)
    remove_empty_responses(conversation.messages)
    remove_consecutive_requests(conversation.messages)
    # Invalidate agent cache - HTTP client may be in bad state after abort
    raise
```

**Implication:** Any partial response is deliberately discarded to prevent corrupted state.

### Pattern 2: Layer Separation via Callback
**Location:** `ui/app.py:348-364` → `core/agents/main.py:388`

Core streams to UI without knowing UI implementation:
```python
# In core layer (main.py:388):
await ac.stream_model_request_node(
    node, agent_run.ctx, state_manager,
    self.streaming_callback,  # ← UI-provided callback
    ...
)

# In UI layer (app.py:354):
async def streaming_callback(self, chunk: str) -> None:
    self.current_stream_text += chunk  # ← Accumulated in UI only
```

**Implication:** Core cannot retrieve accumulated text during abort due to layer boundaries.

### Pattern 3: UI-Level Stream Cancellation
**Location:** `ui/app.py:399-407`

```python
def action_cancel_stream(self) -> None:
    if self._current_request_task is not None:
        self._streaming_cancelled = True
        self._stream_buffer.clear()
        self.current_stream_text = ""  # ← CLEARS before cancel
        self.chat_container.cancel_stream()
        self._current_request_task.cancel()
```

**Implication:** Clearing happens immediately, before core layer can capture partial state.

---

## Knowledge Gaps

1. **Recovery Mechanism**: How should partial responses be marked in the message history? (e.g., `[PARTIAL]` prefix, special part type)

2. **User Experience**: Should partial responses be restored on session load? If so, how should they be displayed?

3. **Resumption**: Should partial responses be included in future agent context? They may confuse the LLM if not properly marked.

4. **Session Format**: The current session JSON format stores `messages` as a list. Adding partial responses requires format changes or new fields.

---

## Potential Solutions

### Option A: Add Partial Response Field to RuntimeState
```python
# In core/types/state_structures.py:
@dataclass(slots=True)
class RuntimeState:
    # ... existing fields ...
    partial_response_text: str = ""  # Store last accumulated stream
    partial_response_timestamp: str = ""
```

**Pros:**
- Clean separation: partial responses tracked separately from message history
- Easy to persist in `save_session()`

**Cons:**
- Requires modifying multiple files (state, UI callback, abort handler)
- Needs UX decision on display

### Option B: Leverage Existing Debug Accumulator
```python
# In core/agents/agent_components/streaming.py:350:
state_manager.session._debug_raw_stream_accum += delta_text

# In abort handler (main.py:438), create partial message:
if state_manager.session._debug_raw_stream_accum:
    partial_text = state_manager.session._debug_raw_stream_accum
    # Create ModelResponse with partial content
    # Add to conversation.messages
```

**Pros:**
- Minimal code changes (accumulator already exists)
- No new state fields

**Cons:**
- `_debug_raw_stream_accum` is intended for debugging, not user content
- Not currently persisted in `save_session()`

### Option C: Bi-Directional Callback Pattern
Add a new callback parameter to retrieve accumulated text from UI:
```python
# Core layer calls UI with both directions:
await streaming_callback(chunk)  # send chunk
partial_text = await get_partial_callback()  # retrieve on abort
```

**Pros:**
- Preserves layer separation (core doesn't access UI directly)
- Flexible for future use cases

**Cons:**
- More complex API surface
- Requires async coordination during abort

---

## References

- **Streaming Callback:** `src/tunacode/ui/app.py:348-364`
- **Cancellation Action:** `src/tunacode/ui/app.py:399-407`
- **Abort Handler:** `src/tunacode/core/agents/main.py:438-465`
- **Message Persistence:** `src/tunacode/core/agents/main.py:206-214`
- **Session Save:** `src/tunacode/core/state.py:246-276`
- **Debug Accumulator:** `src/tunacode/core/agents/agent_components/streaming.py:350`
- **Chat Container:** `src/tunacode/ui/widgets/chat.py:80-150`

**GitHub Permalink:** [8247e85e](https://github.com/alchemiststudiosDOTai/tunacode/blob/8247e85e85f4cead37426260bc819b932472a90a/src/tunacode/core/agents/main.py)
