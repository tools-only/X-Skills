# Research – insert_before_stream() Tool Panel Investigation

**Date:** 2026-01-27
**Owner:** claude-opus-4-5
**Phase:** Research

## Goal

Investigate the `insert_before_stream()` implementation to understand why tool panels might not be appearing correctly during streaming, despite the reported fix.

## Findings

### Core Implementation

**Location:** `/home/tuna/tunacode/src/tunacode/ui/widgets/chat.py:244-261`

```python
def insert_before_stream(self, renderable: RenderableType) -> None:
    """Insert a message before the current streaming widget."""
    widget = Static(renderable)
    widget.add_class("chat-message")

    if self._current_stream is not None:
        # Insert before the streaming widget
        self.mount(widget, before=self._current_stream)
    else:
        # No active stream, just append
        self.mount(widget)

    if self._auto_scroll:
        self.scroll_end(animate=False)
```

**Usage:** `/home/tuna/tunacode/src/tunacode/ui/app.py:283`

```python
# Insert before streaming widget so tool panels appear inline
# without interrupting the stream
self.chat.insert_before_stream(panel)
```

### Flow Analysis

```
Happy Path:
1. start_stream() → _current_stream = MessageWidget
2. Tokens stream → stream(chunk) → append()
3. Tool completes → ToolResultDisplay → insert_before_stream()
   → mount(widget, before=_current_stream) ✓ WORKS
4. end_stream() → finalize() → _current_stream = None
```

### Identified Issues

| Issue | Severity | Description |
|-------|----------|-------------|
| **Cancel + Pending Tools** | HIGH | Tools complete after cancel, `_current_stream` is None, panels append to end |
| **Late Tool Panel** | MEDIUM | Tools complete after `end_stream()`, panels append to bottom |
| **Orphaned Stream** | MEDIUM | Double `start_stream()` silently finalizes first stream |
| **No Throttling** | LOW | Old system had 100ms throttle, new system has none |

### HIGH: Cancel + Pending Tools Race Condition

**Location:** `/home/tuna/tunacode/src/tunacode/ui/app.py:205-207`

```python
except asyncio.CancelledError:
    self.chat.cancel_stream()  # _current_stream = None
    self.notify("Cancelled")
```

**Problem Sequence:**
1. User cancels mid-tool-execution
2. `cancel_stream()` immediately sets `_current_stream = None`
3. Tools are **still running in background** (no cancellation signal)
4. Tool completes → posts `ToolResultDisplay` message
5. `on_tool_result_display()` → `insert_before_stream(panel)`
6. `_current_stream is None` → **panel appends to end instead of inline**

**Visual Result:** Tool panels appear AFTER the "Cancelled" notification, confusing users who expect cancellation to stop everything.

**Evidence:** `cancel_stream()` at chat.py:237-242 does not cancel pending tool executions:
```python
def cancel_stream(self) -> None:
    if self._current_stream is not None:
        self._current_stream.remove()
        self._current_stream = None
    self.remove_class("streaming")
    # No tool cancellation logic
```

### MEDIUM: Late Tool Panel After end_stream()

**Location:** `/home/tuna/tunacode/src/tunacode/ui/app.py:238`

**Problem Sequence:**
1. Model streams rapidly, finishes response
2. `end_stream(panel)` called at line 238
3. `_current_stream = None` (chat.py:231)
4. Slow tool (network delay) completes AFTER this
5. `insert_before_stream()` → `_current_stream is None` → appends to end

**Visual Result:** Tool panel appears BELOW the finalized agent response instead of inline where it belongs chronologically.

### Key Observation

The `insert_before_stream()` design assumes **all tools complete BEFORE the stream ends**. This assumption breaks when:

- Tools have network latency
- User cancels mid-execution
- Parallel tools complete at different times

### Tool Result Lifecycle (Confirmed)

**Flow:**
`core/agents/agent_components/orchestrator/orchestrator.py::_emit_tool_returns()`
→ `tool_result_callback` (built in `ui/repl_support.py::build_tool_result_callback`)
→ `app.post_message(ToolResultDisplay)`
→ `TextualReplApp.on_tool_result_display()`
→ `chat.insert_before_stream(panel)`.

Tool result panels are emitted when **tool-return parts** appear in a node’s request, not when the tool starts. This means delivery timing depends on when the tool return parts are surfaced to the UI loop.

### Textual Mount Ordering (Confirmed)

`Widget.mount(before=...)` resolves a concrete parent + index and calls `App._register(...)`, which inserts into `parent._nodes` via `_insert(before, child)` (Textual source in `.venv/.../textual/widget.py` and `.venv/.../textual/app.py`). Ordering is deterministic **if** the `before` widget is still attached to the same parent; if the streaming widget has been removed, insertion falls back to append in our code.

### Message Queue Timing (Confirmed)

`MessagePump.post_message()` **queues** the message into an asyncio queue and returns immediately; handling happens later in `_process_messages_loop()` (`.venv/.../textual/message_pump.py`). Messages are not coalesced by default (`Message.can_replace()` returns False), so `ToolResultDisplay` is processed asynchronously relative to stream completion.

## Related Files

- `src/tunacode/ui/widgets/chat.py` → ChatContainer, MessageWidget, insert_before_stream
- `src/tunacode/ui/app.py:271-283` → on_tool_result_display handler
- `src/tunacode/ui/app.py:205-207` → Cancel handling
- `src/tunacode/core/agents/main.py:438-465` → Abort cleanup (message history only, no UI sync)
- `memory-bank/research/2026-01-27_17-45-00_chatcontainer_refactor.md` → Previous refactor research

## Knowledge Gaps (Resolved)

1. **Tool execution lifecycle** - Tool results are emitted from `_emit_tool_returns()` during node processing; `build_tool_result_callback()` posts `ToolResultDisplay`, handled by `TextualReplApp.on_tool_result_display()`.
2. **Textual mount ordering** - `mount(before=widget)` inserts into the parent’s node list at a specific index; ordering is deterministic while the reference widget remains mounted.
3. **Message queue timing** - `post_message()` queues to the message pump; delivery is async and can occur after stream end/cancel.

## References

- `chat.py:244-261` - insert_before_stream implementation
- `app.py:205-207` - Cancel handling
- `app.py:271-283` - Tool result display handler
- Previous research: `memory-bank/research/2026-01-27_17-45-00_chatcontainer_refactor.md`
