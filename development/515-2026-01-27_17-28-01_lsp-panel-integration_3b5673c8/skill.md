# Research – LSP Panel Integration (Dependency-Compliant)

**Date:** 2026-01-27 17:28:01
**Owner:** agent
**Phase:** Research

## Goal

Research how to add LSP status/diagnostics back to a UI panel after the refactor, while strictly following the dependency direction rule: `ui → core → tools → utils/types`.

## Context

LSP was recently refactored (ticket `tun-c2f8`):
- **Before:** LSP lived in `src/tunacode/lsp/` with status indicators in UI
- **After:** LSP moved to `src/tunacode/tools/lsp/` (tools-only concern)
- **Removed:** LSP status indicator from `ResourceBar`

The user wants to restore LSP status in the UI without violating dependency boundaries.

## Findings

### Current LSP Architecture

**Location:** `src/tunacode/tools/lsp/`

| File | Purpose |
|------|---------|
| `__init__.py` | Public API: `get_diagnostics()`, `format_diagnostics()` |
| `client.py` | `LSPClient` class, `Diagnostic` dataclass |
| `servers.py` | Language server command mappings |
| `diagnostics.py` | `maybe_prepend_lsp_diagnostics()`, `is_lsp_enabled()` |

**Current Integration:**
- `write_file.py` and `update_file.py` call `maybe_prepend_lsp_diagnostics()`
- Diagnostics are embedded in tool result as `<file_diagnostics>` block
- UI renderer (`ui/renderers/tools/diagnostics.py`) parses and displays inline

### Existing Panel Components

**Status Displays:**
- `ResourceBar` (`ui/widgets/resource_bar.py`) - tokens, costs, model
- `StatusBar` (`ui/widgets/status_bar.py`) - location, edited files, last action

**LSP was previously in ResourceBar** but was removed during refactor.

### Dependency-Compliant Integration Patterns

The codebase uses 4 patterns for UI ← Core/Tools data flow:

#### Pattern 1: Callback Injection (Streaming, Tool Results)
```python
# UI builds callback
def build_tool_result_callback(app: AppForCallbacks) -> ToolResultCallback:
    def _callback(tool_name, status, args, result, duration_ms):
        app.post_message(ToolResultDisplay(...))
    return _callback

# UI passes to core
await process_request(..., tool_result_callback=callback)

# Core invokes callback (never imports UI)
callback(tool_name, status, args, result, duration_ms)
```

#### Pattern 2: Shared State via StateManager (ResourceBar)
```python
# UI holds reference
self.state_manager: StateManager = state_manager

# Core/Tools write to state
state_manager.session.some_value = new_value

# UI reads from state
def _update_resource_bar(self) -> None:
    session = self.state_manager.session
    self.resource_bar.update_stats(tokens=session.total_tokens)
```

#### Pattern 3: Protocol-Based Decoupling
```python
# UI defines protocol (repl_support.py)
class StatusBarLike(Protocol):
    def update_last_action(self, tool_name: str) -> None: ...

class AppForCallbacks(Protocol):
    status_bar: StatusBarLike
    def post_message(self, message: ToolResultDisplay) -> bool: ...

# Core depends on protocol, not concrete class
```

#### Pattern 4: Textual Message Bus
```python
# Callback posts message
app.post_message(ToolResultDisplay(...))

# Handler renders
def on_tool_result_display(self, message: ToolResultDisplay) -> None:
    panel = tool_panel_smart(...)
    self.chat.add_message(panel)
```

## Recommended Integration Approach

### Option A: StateManager Pattern (Recommended)

**Cleanest approach - follows existing ResourceBar pattern.**

1. **Add LSP state to Session** (`core/state.py` or `core/types.py`):
   ```python
   @dataclass
   class LSPStatus:
       enabled: bool = False
       last_check: datetime | None = None
       diagnostics_count: int = 0
       server_running: bool = False
   ```

2. **Tools write to state** (`tools/lsp/diagnostics.py`):
   ```python
   def maybe_prepend_lsp_diagnostics(..., state_manager=None):
       # ... existing logic ...
       if state_manager:
           state_manager.session.lsp_status.diagnostics_count = len(diagnostics)
   ```

3. **UI reads from state** (`ui/app.py`):
   ```python
   def _update_resource_bar(self) -> None:
       session = self.state_manager.session
       self.resource_bar.update_stats(
           tokens=session.total_tokens,
           lsp_enabled=session.lsp_status.enabled,
           lsp_diagnostics=session.lsp_status.diagnostics_count,
       )
   ```

**Dependency flow:** `ui → core.state` (clean)

### Option B: Callback Pattern

**More explicit, but adds threading through layers.**

1. **Add LSP callback type** (`core/types.py`):
   ```python
   LSPStatusCallback = Callable[[bool, int], None]  # (enabled, diagnostic_count)
   ```

2. **UI builds callback** (`ui/repl_support.py`):
   ```python
   def build_lsp_status_callback(app: AppForCallbacks) -> LSPStatusCallback:
       def _callback(enabled: bool, count: int) -> None:
           app.resource_bar.update_lsp_status(enabled, count)
       return _callback
   ```

3. **Thread through process_request** → agent → tool execution

**Dependency flow:** `ui builds callback → core → tools invoke callback` (clean but verbose)

### Option C: Diagnostic Extraction from Tool Results

**Minimal change - extract from existing data.**

1. **UI extracts from tool result** (`ui/app.py`):
   ```python
   def on_tool_result_display(self, message: ToolResultDisplay) -> None:
       # Extract diagnostics from result
       diagnostics = extract_diagnostics_from_result(message.result)
       if diagnostics:
           self.resource_bar.update_lsp_indicator(len(diagnostics.items))
   ```

**Dependency flow:** `ui → ui.renderers.tools.diagnostics` (internal, clean)

## Key Patterns / Solutions Found

| Pattern | Use Case | Example Location |
|---------|----------|------------------|
| Callback injection | Real-time events | `repl_support.py:147-174` |
| StateManager shared state | Persistent status | `app.py:389-400` |
| Protocol decoupling | Type-safe contracts | `repl_support.py:95-107` |
| Textual messages | Widget communication | `widgets/messages.py` |

## Knowledge Gaps (Resolved)

1. ~~**StateManager threading** - Does `StateManager` get passed to tool functions currently?~~
   **ANSWER: NO.** Tools don't receive `state_manager`. This rules out Option A without significant refactoring.

2. **LSP lifecycle** - Updates per-file (when `write_file`/`update_file` complete)

3. **ResourceBar capacity** - TBD (visual design needed)

4. **Performance** - Not a concern with Option C (only updates on tool completion)

## Recommendation: Option C (Extraction)

**Why:** Infrastructure already exists, no architectural changes needed.

**Existing functions in `ui/renderers/tools/diagnostics.py`:**
- `parse_diagnostics_block(result)` → `DiagnosticsData` with counts
- `extract_diagnostics_from_result(result)` → separates diagnostics from result
- `DiagnosticsData.error_count`, `.warning_count`, `.info_count` already computed

## Implementation Checklist (Option C - Recommended)

**Minimal changes, uses existing infrastructure:**

1. **ResourceBar** (`ui/widgets/resource_bar.py`):
   - [ ] Add `_lsp_errors: int = 0`, `_lsp_warnings: int = 0` state
   - [ ] Add `update_lsp_status(errors: int, warnings: int)` method
   - [ ] Update `_refresh_display()` to show LSP indicator

2. **App handler** (`ui/app.py:on_tool_result_display()`):
   - [ ] Import `parse_diagnostics_block` from `ui/renderers/tools/diagnostics`
   - [ ] Extract diagnostics: `data = parse_diagnostics_block(message.result)`
   - [ ] If data: `self.resource_bar.update_lsp_status(data.error_count, data.warning_count)`

3. **Optional - Session tracking**:
   - [ ] Add `_session_lsp_errors: int` to track cumulative
   - [ ] Reset on new session

**Dependency flow:** `ui/app.py` → `ui/renderers/tools/diagnostics.py` (same layer, clean)

## References

**Files:**
- `src/tunacode/tools/lsp/` - Current LSP implementation
- `src/tunacode/ui/widgets/resource_bar.py` - Target panel
- `src/tunacode/ui/repl_support.py` - Callback patterns
- `src/tunacode/core/state.py` - StateManager
- `src/tunacode/ui/renderers/tools/diagnostics.py` - Existing diagnostic parsing

**Tickets:**
- `.tickets/tun-c2f8.md` - LSP refactor epic (completed)

**Documentation:**
- `docs/lsp-diagnostics.md` - User-facing LSP docs
- `docs/codebase-map/modules/lsp.md` - LSP architecture
