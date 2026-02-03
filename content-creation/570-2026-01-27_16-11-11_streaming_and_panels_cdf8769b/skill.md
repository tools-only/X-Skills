# Research – Streaming Response and Panel Display in TunaCode TUI

**Date:** 2026-01-27 16:11:11
**Owner:** Claude
**Phase:** Research
**Git Commit:** 601a7c1c
**Last Updated:** 2026-01-27 16:11:11

## Goal

Explain how streaming responses work in the TunaCode TUI and what panels are shown before the final answer. This includes understanding the dual-display architecture (streaming overlay vs. permanent RichLog), the callback chain from model to display, panel rendering system, and the visual state management during active generation.

---

## Executive Summary

TunaCode implements a **dual-display system** for model responses:

1. **Real-time streaming panel** (`#streaming-output`): A fixed-position overlay that updates incrementally as tokens arrive, positioned between the viewport and editor
2. **Permanent RichLog panels**: Finalized responses and tool results written to the scrollable history

The key insight is that **streaming content is NOT written directly to RichLog**. Instead, it accumulates in a separate Static widget that appears only during active generation. This prevents scroll jitter and keeps the active thought process focused. When streaming completes, the temporary panel is cleared and a permanent panel is written to RichLog with metadata (token count, duration, timestamp).

---

## Findings

### Relevant Files & Why They Matter

| File | Purpose |
|------|---------|
| `src/tunacode/core/agents/agent_components/streaming.py` | Core streaming instrumentation, token-level streaming from LLM provider, callback invocation |
| `src/tunacode/ui/app.py` (lines 328-344, 346-351) | Main TUI app: `streaming_callback()`, `_update_streaming_panel()`, streaming_output widget management |
| `src/tunacode/ui/renderers/agent_response.py` | Panel renderers: `render_agent_streaming()` (active) vs `render_agent_response()` (final) |
| `src/tunacode/ui/renderers/panels.py` | Generic panel rendering system, PanelType enum, color schemes |
| `src/tunacode/types/callbacks.py` | StreamingCallback type definition (contract for streaming) |
| `src/tunacode/ui/styles/layout.tcss` | CSS for #streaming-output (fixed dock, .active class for visibility) |
| `src/tunacode/ui/styles/theme-nextstep.tcss` | NeXTSTEP theme styling for streaming panels (inverted bevel effect) |

### Architectural Diagram

```
User Input (Editor)
    ↓
_request_worker() queue
    ↓
_process_request() [app.py:181-245]
    ↓
process_request() [main.py:536-556]
    ↓
RequestOrchestrator.run() [main.py:216-234]
    ↓
async with agent.iter() [main.py:350]
    ↓
stream_model_request_node() [streaming.py:122]
    ↓
async for event in request_stream [streaming.py:209]
    ↓
PartDeltaEvent → TextPartDelta.content_delta [streaming.py:283-351]
    ↓
await streaming_callback(delta_text) [streaming.py:351]
    ↓
┌─────────────────────────────────────────────────┐
│ DUAL DISPLAY SYSTEM                             │
├─────────────────────────────────────────────────┤
│ 1. STREAMING PHASE (active):                    │
│    streaming_callback() [app.py:328]            │
│    → current_stream_text += chunk               │
│    → throttle check (100ms)                     │
│    → _update_streaming_panel() [app.py:346]     │
│    → render_agent_streaming() [agent_response.py:86]
│    → streaming_output.update(panel)             │
│    → [Shows in #streaming-output overlay]       │
├─────────────────────────────────────────────────┤
│ 2. FINAL PHASE (completion):                    │
│    streaming_output.update("")                  │
│    → render_agent_response() [agent_response.py:152]
│    → rich_log.write(panel, expand=True)         │
│    → [Permanent panel in RichLog]               │
└─────────────────────────────────────────────────┘
```

---

## Key Patterns / Solutions Found

### 1. Streaming Callback Architecture

**Contract** (`src/tunacode/types/callbacks.py:37`):
```python
StreamingCallback: TypeAlias = Callable[[str], Awaitable[None]]
```

**Preconditions**: chunk is ordered text delta.
**Postconditions**: enqueue or render the chunk without raising.

The callback is passed through the agent stack:
1. Created in `ui/app.py` as a bound method
2. Passed to `core/agents/main.py` via `process_request()`
3. Passed to `stream_model_request_node()` in `streaming.py`
4. Invoked per token chunk: `await streaming_callback(delta_text)`

### 2. Throttled Display Updates

To reduce visual churn and prevent excessive markdown re-renders:

```python
STREAM_THROTTLE_MS: float = 100.0  # [app.py:54]

async def streaming_callback(self, chunk: str) -> None:
    self.current_stream_text += chunk  # Always accumulate immediately

    now = time.monotonic()
    elapsed_ms = (now - self._last_display_update) * 1000

    if elapsed_ms >= STREAM_THROTTLE_MS:
        self._last_display_update = now
        self._update_streaming_panel(now)
        self.streaming_output.add_class("active")
        self.rich_log.scroll_end()
```

### 3. Streaming Panel vs. Final Panel Differences

| Aspect | Streaming Panel | Final Panel |
|--------|----------------|-------------|
| Renderer | `render_agent_streaming()` [agent_response.py:86] | `render_agent_response()` [agent_response.py:152] |
| Widget | `#streaming-output` (Static) | `rich_log` (RichLog) |
| Title indicator | `[...]` shows active | Timestamp shows completion |
| Border color | `primary` (pink/cyan) | `accent` |
| Metadata | Elapsed time | Tokens, tokens/sec, duration |
| Visibility | Temporary (`.active` class) | Permanent (written to RichLog) |
| Position | Fixed dock between viewport/editor | Scrollable in viewport |

### 4. CSS Class Management for Visual States

| Class | Widget | Effect | Added By | Removed By |
|-------|--------|--------|----------|------------|
| `streaming` | `#viewport` | Accent border during generation | `_process_request()` line 185 | `finally` line 218 |
| `paused` | `#viewport` | Double border during pause | `pause_streaming()` line 361 | `resume_streaming()` line 366 |
| `active` | `#streaming-output` | Shows the streaming panel | `streaming_callback()` line 343 | `finally` line 221 |

### 5. Tool Panel Routing

Tool results use a smart routing system (`src/tunacode/ui/renderers/panels.py:479-534`):

```python
def tool_panel_smart(name, status, args, result, ...):
    # Only completed tools with results get custom renderers
    if status == "completed" and result:
        renderer_map = {
            "list_dir": render_list_dir,
            "grep": render_grep,
            "glob": render_glob,
            "read_file": render_read_file,
            "update_file": render_update_file,
            "bash": render_bash,
            "web_fetch": render_web_fetch,
        }
        renderer = renderer_map.get(name.lower())
        if renderer:
            return renderer(args, result, duration_ms, max_line_width)

    # Fallback to generic panel
    return tool_panel(name, status, args, result, ...)
```

Tool panels are written directly to RichLog via `ToolResultDisplay` message events, NOT to the streaming overlay.

### 6. Pause/Resume Buffering

When user presses `Ctrl+P` to pause streaming:

```python
self._streaming_paused: bool = False
self._stream_buffer: list[str] = []  # [app.py:80]

async def streaming_callback(self, chunk: str) -> None:
    if self._streaming_paused:
        self._stream_buffer.append(chunk)  # Buffer instead of accumulate
        return
    # ... normal accumulation and display
```

On resume, buffer is flushed:
```python
def resume_streaming(self) -> None:
    if self._stream_buffer:
        buffered_text = "".join(self._stream_buffer)
        self.current_stream_text += buffered_text  # Flush buffer
        self._stream_buffer.clear()

    # Force immediate display update on resume
    self._update_streaming_panel(now)
```

### 7. NeXTSTEP Panel Layout

All panels follow a 4-zone NeXTSTEP layout:

```
┌─ Header: tool_name [status_suffix] ─────────────┐
│ Parameters: selection context (args)             │
├──────────────────────────────────────────────────┤
│ Viewport: primary content                        │
│ (markdown, syntax-highlighted code, etc.)        │
├──────────────────────────────────────────────────┤
│ Status: metrics, truncation info, timestamp      │
└──────────────────────────────────────────────────┘
```

See `docs/ui/nextstep_panels.md` for detailed guidelines.

---

## Behavioral States

### User Request Lifecycle

1. **Idle State**: No request in progress
   - `#streaming-output` hidden (`display: none`)
   - Editor focused, awaiting input

2. **Loading State**: Request queued, awaiting first token
   - `LoadingIndicator` shows spinner
   - Viewport has no classes yet

3. **Streaming State**: Tokens arriving, display updating
   - `#viewport` has `.streaming` class (accent border)
   - `#streaming-output.active` shows streaming panel
   - Panel updates every 100ms (throttled)
   - `[...]` indicator in title shows active state

4. **Paused State**: User pressed Ctrl+P
   - `#viewport` has `.paused` class (double border)
   - Chunks buffered, display not updated
   - Resume flushes buffer and forces immediate update

5. **Completion State**: Streaming finished
   - `#streaming-output` cleared and hidden
   - Final panel written to RichLog with metadata
   - All CSS classes removed

6. **Cancellation State**: User pressed Escape
   - `_streaming_cancelled = True`
   - Streaming panel cleared
   - No final panel written to RichLog

---

## Design Decision: Why Separate Streaming Panel?

From `docs/ui/layout.md` (lines 44-49):

> **Streaming Output** - A specialized dock that appears between the history and the editor. Prevents the "jumpy" scrolling effect that occurs when streaming text directly into a long log. It keeps the active thought process static and focused.

**Technical rationale**:
- RichLog is a scrollable container with potentially hundreds of panels
- Updating content inside RichLog causes scroll position jitter
- By separating active streaming into a fixed sibling widget, we:
  - Keep the stream visible without scrolling
  - Maintain smooth RichLog performance
  - Provide clear visual separation between "in progress" and "done"

---

## Knowledge Gaps

- [ ] What happens when streaming callback raises an exception?
- [ ] How does the system handle malformed markdown during streaming?
- [ ] Are there any race conditions between throttled updates and completion?
- [ ] How does the panel width calculation handle terminal resize events?

---

## References

### Implementation Files

- [src/tunacode/core/agents/agent_components/streaming.py](https://github.com/alchemiststudiosDOTai/tunacode/blob/601a7c1c/src/tunacode/core/agents/agent_components/streaming.py) - Core streaming orchestration
- [src/tunacode/ui/app.py:328-351](https://github.com/alchemiststudiosDOTai/tunacode/blob/601a7c1c/src/tunacode/ui/app.py#L328-L351) - Streaming callback and panel update
- [src/tunacode/ui/renderers/agent_response.py](https://github.com/alchemiststudiosDOTai/tunacode/blob/601a7c1c/src/tunacode/ui/renderers/agent_response.py) - Panel renderers
- [src/tunacode/ui/renderers/panels.py](https://github.com/alchemiststudiosDOTai/tunacode/blob/601a7c1c/src/tunacode/ui/renderers/panels.py) - Generic panel system
- [src/tunacode/types/callbacks.py](https://github.com/alchemiststudiosDOTai/tunacode/blob/601a7c1c/src/tunacode/types/callbacks.py) - Callback type definitions

### Styles

- [src/tunacode/ui/styles/layout.tcss:68-81](https://github.com/alchemiststudiosDOTai/tunacode/blob/601a7c1c/src/tunacode/ui/styles/layout.tcss#L68-L81) - Streaming output CSS
- [src/tunacode/ui/styles/theme-nextstep.tcss:105-120](https://github.com/alchemiststudiosDOTai/tunacode/blob/601a7c1c/src/tunacode/ui/styles/theme-nextstep.tcss#L105-L120) - NeXTSTEP panel styling

### Documentation

- [docs/ui/layout.md](https://github.com/alchemiststudiosDOTai/tunacode/blob/601a7c1c/docs/ui/layout.md) - Main layout specification
- [docs/ui/agent_response_panel.md](https://github.com/alchemiststudiosDOTai/tunacode/blob/601a7c1c/docs/ui/agent_response_panel.md) - Agent response panel specification
- [docs/ui/nextstep_panels.md](https://github.com/alchemiststudiosDOTai/tunacode/blob/601a7c1c/docs/ui/nextstep_panels.md) - NeXTSTEP panel guidelines
- [docs/ui/tool_renderers.md](https://github.com/alchemiststudiosDOTai/tunacode/blob/601a7c1c/docs/ui/tool_renderers.md) - Tool renderer patterns

### Debug History

- [.claude/debug_history/2026-01-21_stream-hang-timeout.md](https://github.com/alchemiststudiosDOTai/tunacode/blob/601a7c1c/.claude/debug_history/2026-01-21_stream-hang-timeout.md) - Stream hang investigation
- [.claude/debug_history/2026-01-21_resume-hang-investigation.md](https://github.com/alchemiststudiosDOTai/tunacode/blob/601a7c1c/.claude/debug_history/2026-01-21_resume-hang-investigation.md) - Resume hang investigation
