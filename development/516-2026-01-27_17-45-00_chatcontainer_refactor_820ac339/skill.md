# Research – ChatContainer + MessageWidget Refactoring

**Date:** 2026-01-27 17:45:00
**Owner:** Claude
**Phase:** Research
**Git Commit:** 601a7c1c
**Last Updated:** 2026-01-27 17:45:00

## Goal

Validate the proposed refactoring to replace the current RichLog + streaming_output dual-display system with a unified ChatContainer + MessageWidget architecture. Identify all touch points, dependencies, and risks.

---

## Executive Summary

The proposed refactoring is **architecturally sound** but involves **33 touch points across 4 files**. The current system uses:
- `RichLog` as the permanent chat history (write-only, append-based)
- `streaming_output` (Static widget) as an ephemeral overlay for live streaming

The proposed `ChatContainer + MessageWidget` unifies these into a single scrollable container where messages can be updated in-place during streaming, then finalized. This eliminates the dual-display complexity but requires careful migration of all write sites.

---

## Findings

### Current Architecture

```
┌─────────────────────────────────────────────────────────┐
│ ResourceBar (fixed)                                      │
├─────────────────────────────────────────────────────────┤
│ #viewport Container                                      │
│   ├── RichLog (scrollable, permanent history)           │
│   └── LoadingIndicator (overlay)                        │
├─────────────────────────────────────────────────────────┤
│ #streaming-output Static (ephemeral, docked)            │
├─────────────────────────────────────────────────────────┤
│ Editor (fixed)                                           │
│ AutoComplete widgets (overlay)                          │
├─────────────────────────────────────────────────────────┤
│ StatusBar (fixed)                                        │
└─────────────────────────────────────────────────────────┘
```

**Key Files:**
| File | Role |
|------|------|
| `src/tunacode/ui/app.py` | Main TUI app, owns RichLog and streaming_output |
| `src/tunacode/ui/renderers/agent_response.py` | Renders streaming vs finalized panels |
| `src/tunacode/ui/renderers/panels.py` | Tool panel rendering system |
| `src/tunacode/ui/commands/__init__.py` | Command output to RichLog |
| `src/tunacode/ui/welcome.py` | Welcome screen to RichLog |

### Existing Widget Inventory

**Location:** `src/tunacode/ui/widgets/`

| Widget | Base Class | Purpose |
|--------|------------|---------|
| `Editor` | `Input` | Main text input with paste buffer |
| `ResourceBar` | `Static` | Model/token/context display |
| `StatusBar` | `Horizontal` | Action status display |
| `CommandAutoComplete` | `AutoComplete` | Slash command completion |
| `FileAutoComplete` | `AutoComplete` | File path completion |

**No existing chat/message widgets** - the codebase uses Textual's `RichLog` directly.

### Touch Points Requiring Changes

#### app.py (24 locations)

**Widget Declaration & Instantiation:**
```python
# Line 90 - Type annotation
self.rich_log: RichLog

# Line 98 - Constructor
self.rich_log = RichLog(wrap=True, markup=True, highlight=True, auto_scroll=True)

# Line 106 - Composition
yield self.rich_log
```

**Write Operations (21 sites):**
| Line | Content Type | Method |
|------|--------------|--------|
| 166 | Debug logs | `write(renderable)` |
| 182 | Error panel | `write(error_renderable)` |
| 218 | Error panel | `write(error_renderable)` |
| 242 | Blank line | `write("")` |
| 243 | Agent response | `write(panel, expand=True)` |
| 264 | Blank line | `write("")` |
| 270 | User message | `write(user_block)` |
| 282 | Tool panel | `write(panel, expand=True)` |
| 301 | System notice | `write(notice_text)` |
| 328 | Replayed user | `write(user_block)` |
| 330 | Replayed label | `write(Text("agent:"))` |
| 331 | Replayed agent | `write(Markdown(content))` |
| 404 | Shell output | `write(renderable)` |

**Width Access (3 sites):**
```python
# Line 265 - User message width
render_width = max(1, self.rich_log.size.width - 2)

# Lines 287, 289 - Tool panel width
self.rich_log.content_region.width
self.rich_log.size.width
```

**Scroll Control (1 site):**
```python
# Line 349 - During streaming
self.rich_log.scroll_end()
```

#### commands/__init__.py (11 locations)

| Line | Content | Purpose |
|------|---------|---------|
| 69 | Table | Help command |
| 79 | `clear()` | Clear command |
| 132 | Text | Debug status |
| 163 | Text | Model config error |
| 376 | `clear()` | Resume clear |
| 385 | Text | Resume loaded message |
| 411-414 | Text | Version info (4 lines) |
| 454 | Text | Update success |
| 458 | Text | Update error |
| 460 | Text | Update exception |

#### welcome.py (3 locations)

| Line | Content | Purpose |
|------|---------|---------|
| 32 | Text | Logo load error |
| 34 | Logo | ASCII art logo |
| 61 | Text | Welcome message |

### Streaming Callback Flow

```
LLM Provider → pydantic-ai events
    ↓
streaming.py:351 → await streaming_callback(delta_text)
    ↓
app.py:333 → streaming_callback(chunk)
    ↓
app.py:339 → self.current_stream_text += chunk (immediate accumulation)
    ↓
app.py:345-350 → throttle check (100ms)
    ↓
app.py:356 → render_agent_streaming() → streaming_output.update(panel)
    ↓
[On completion]
    ↓
app.py:236-243 → render_agent_response() → rich_log.write(panel, expand=True)
    ↓
app.py:225-226 → streaming_output.update("") + remove_class("active")
```

**Key Insight:** Streaming uses `Static.update()` (in-place replacement), not `RichLog.write()` (append). The proposed `MessageWidget.append()` pattern matches this better.

### API Surface Required for ChatContainer

```python
class ChatContainer(ScrollableContainer):
    """Must provide these APIs for compatibility."""

    # Core write API
    def write(self, renderable: RenderableType, *, expand: bool = False) -> None: ...

    # Dimension APIs
    @property
    def size(self) -> Size: ...  # Total widget dimensions

    @property
    def content_region(self) -> Region: ...  # Available content area

    # Control APIs
    def clear(self) -> None: ...
    def scroll_end(self, *, animate: bool = True) -> None: ...

    # New streaming APIs
    def start_assistant_message(self) -> MessageWidget: ...
    def stream_chunk(self, chunk: str) -> None: ...
```

### Width Calculation Anti-Pattern

**Critical:** `expand=True` does NOT mean terminal width!

```python
# WRONG assumption
panel_width = terminal_width  # NO!

# CORRECT
panel_width = scrollable_content_region.width  # Excludes padding + scrollbar (4-8 chars narrower)
```

The proposed refactoring must preserve this behavior. `MessageWidget` should use its parent container's `content_region.width`, not `app.size.width`.

---

## Proposed Architecture

```
┌─────────────────────────────────────────────────────────┐
│ ResourceBar (fixed)                                      │
├─────────────────────────────────────────────────────────┤
│ #viewport Container                                      │
│   └── ChatContainer (id="chat-container")               │
│         ├── MessageWidget(role="user")                  │
│         ├── MessageWidget(role="assistant", streaming)  │
│         ├── MessageWidget(role="tool")                  │
│         ├── MessageWidget(role="system")                │
│         └── LoadingIndicator (overlay)                  │
├─────────────────────────────────────────────────────────┤
│ Editor (fixed)                                           │
│ AutoComplete widgets (overlay)                          │
├─────────────────────────────────────────────────────────┤
│ StatusBar (fixed)                                        │
└─────────────────────────────────────────────────────────┘
```

**Key Change:** `#streaming-output` is eliminated. Streaming happens in-place within `MessageWidget`.

### Migration Mapping

| Current Pattern | New Pattern |
|-----------------|-------------|
| `rich_log.write(user_block)` | `chat.add_message("user", content)` |
| `rich_log.write(panel, expand=True)` | `chat.add_message("assistant").set_content(panel)` |
| `streaming_output.update(panel)` | `chat._current_message.append(chunk)` |
| `rich_log.write(tool_panel)` | `chat.add_message("tool").set_content(panel)` |
| `rich_log.write(notice)` | `chat.add_message("system", notice)` |
| `rich_log.clear()` | `chat.clear()` |
| `rich_log.scroll_end()` | `chat.scroll_end()` |

---

## Effort Estimate

| Task | Files Changed | Complexity | Risk |
|------|---------------|------------|------|
| Create `widgets/chat.py` | 1 new | Low | Low |
| Refactor `app.py` compose/streaming | 1 | Medium | Medium |
| Update message handlers in `app.py` | 1 | Low | Low |
| Update `commands/__init__.py` | 1 | Low | Low |
| Update `welcome.py` | 1 | Low | Low |
| Update CSS (`.message-*` classes) | 1 | Low | Low |
| Update session replay | 1 | Low | Low |
| Testing & edge cases | N/A | Medium | Medium |

**Total:** ~1 day focused work (as proposed in the ticket)

---

## Risks & Mitigations

### Risk 1: Width Calculation Regressions

**Risk:** Tool panels render too narrow or too wide after migration.

**Mitigation:**
- Preserve `tool_panel_max_width()` logic in `app.py`
- Pass explicit width to `MessageWidget._render_content()` rather than relying on `expand=True`

### Risk 2: Scroll Behavior Changes

**Risk:** Auto-scroll during streaming breaks or becomes janky.

**Mitigation:**
- `ChatContainer.stream_chunk()` must call `scroll_end()` after each update
- Test with long responses that exceed viewport height

### Risk 3: Performance Degradation

**Risk:** Creating new `MessageWidget` for each message is slower than `RichLog.write()`.

**Mitigation:**
- `MessageWidget` should use `Static.update()` internally (same as current streaming_output)
- Avoid re-mounting widgets - update content in-place
- Maintain 100ms throttle on display updates

### Risk 4: Backward Compatibility

**Risk:** Commands that access `app.rich_log` directly will break.

**Mitigation:**
- Keep `self.rich_log` as an alias to `self.chat` during transition, or
- Update all command files to use new `chat` API

---

## Knowledge Gaps

- [ ] How does `MessageWidget.append()` affect screen reader accessibility?
- [ ] What's the memory impact of keeping many `MessageWidget` instances vs RichLog internal buffer?
- [ ] Does `ScrollableContainer.scroll_end()` have same signature as `RichLog.scroll_end()`?

---

## References

### Implementation Files
- `src/tunacode/ui/app.py:90-424` - Main TUI app with RichLog usage
- `src/tunacode/ui/renderers/agent_response.py:86-223` - Streaming vs finalized renderers
- `src/tunacode/ui/renderers/panels.py:479-534` - Tool panel smart routing
- `src/tunacode/ui/commands/__init__.py` - Command RichLog writes
- `src/tunacode/ui/welcome.py:27-61` - Welcome screen RichLog writes

### Related Research
- `memory-bank/research/2026-01-27_16-11-11_streaming_and_panels.md` - Prior streaming architecture research

### Widget References
- `src/tunacode/ui/widgets/__init__.py` - Current widget exports
- `src/tunacode/ui/widgets/editor.py` - Example of custom widget extending Textual base

### CSS Files
- `src/tunacode/ui/styles/layout.tcss` - RichLog and streaming-output styles
- `src/tunacode/ui/styles/panels.tcss` - Panel styling
- `src/tunacode/ui/styles/theme-nextstep.tcss` - NeXTSTEP theme

---

## Recommendations

1. **Start with the minimal implementation** from the proposal:
   ```python
   class MessageWidget(Static):
       def append(self, chunk: str) -> None:
           self.content += chunk
           self.update(Panel(Markdown(self.content), title=self.role))
   ```

2. **Add width parameter** to `_render_content()`:
   ```python
   def _render_content(self, width: int | None = None) -> None:
       panel = render_agent_streaming(self.content, 0, "", width=width)
       self.update(panel)
   ```

3. **Preserve compatibility layer** initially:
   ```python
   # In app.py
   @property
   def rich_log(self) -> ChatContainer:
       """Backward compat alias."""
       return self.chat
   ```

4. **Update tests** in parallel with implementation to catch regressions early.

5. **Follow NeXTSTEP-ui skill** for any visual changes - the MessageWidget styling must match existing panel aesthetics.
