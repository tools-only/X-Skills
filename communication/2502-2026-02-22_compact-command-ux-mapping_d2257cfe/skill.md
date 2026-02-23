# Research – /compact Command UI/UX Mapping

**Date:** 2026-02-22
**Owner:** claude
**Phase:** Research
**git_commit:** 050ebe24f3918159f2de53544a3eb2e445a94922

## Goal

Map the complete UI/UX flow of the `/compact` command, identify pain points, and catalog the feedback mechanisms available to the user during compaction.

## Findings

### Data Flow: Controller → Callback → UI

```
CompactCommand.execute()
  ├─ controller.set_status_callback(app._update_compaction_status)     # compact.py:53
  ├─ controller.force_compact(...)                                      # compact.py:58
  │    └─ _compact()                                                    # controller.py:210
  │         ├─ _announce_compaction_start() → callback(True)            # controller.py:236
  │         ├─ summarizer.summarize()                                   # controller.py:238
  │         └─ finally: _announce_compaction_end() → callback(False)    # controller.py:264
  ├─ finally: controller.set_status_callback(None)                      # compact.py:64
  ├─ finally: app._update_compaction_status(False)                      # compact.py:65
  └─ finally: app._update_resource_bar()                                # compact.py:66
```

### Relevant Files

| File | Purpose |
|------|---------|
| [`src/tunacode/ui/commands/compact.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/050ebe24/src/tunacode/ui/commands/compact.py) | `/compact` slash command handler |
| [`src/tunacode/core/compaction/controller.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/050ebe24/src/tunacode/core/compaction/controller.py) | Compaction orchestrator, callback registration |
| [`src/tunacode/core/compaction/types.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/050ebe24/src/tunacode/core/compaction/types.py) | Status types: `compacted`, `skipped`, `failed` |
| [`src/tunacode/ui/app.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/050ebe24/src/tunacode/ui/app.py) | `_update_compaction_status` (L531), `_update_resource_bar` (L534), `action_toggle_context_panel` (L430) |
| [`src/tunacode/ui/widgets/resource_bar.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/050ebe24/src/tunacode/ui/widgets/resource_bar.py) | ResourceBar widget, renders "Compacting..." label |
| [`src/tunacode/ui/context_panel.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/050ebe24/src/tunacode/ui/context_panel.py) | Context panel (no compaction status display) |
| [`src/tunacode/ui/styles/layout.tcss`](https://github.com/alchemiststudiosDOTai/tunacode/blob/050ebe24/src/tunacode/ui/styles/layout.tcss) | `ResourceBar.hidden { display: none }` (L33-35) |

### Current User Feedback During Compaction

| Phase | Feedback Type | Location |
|-------|--------------|----------|
| Args validation fail | Toast notification | `compact.py:35` |
| Empty history | Toast notification | `compact.py:48` |
| During compaction | "Compacting..." text in ResourceBar | `resource_bar.py:127-129` |
| Skip/no-op | Toast notification (warning) | `compact.py:78-81` |
| Success | Toast notification with stats | `compact.py:87-92` |
| Error/exception | Nothing visible (silently caught in `finally`) | `compact.py:63-66` |

### Ctrl+E Interaction (The Spaz Bug Context)

The ResourceBar and context panel are **mutually exclusive** via CSS `display: none`:

```python
# app.py:419-428
def _set_context_panel_visibility(self, *, visible: bool) -> None:
    if visible:
        context_rail.remove_class("hidden")      # show panel
        self.resource_bar.add_class("hidden")     # HIDE resource bar
    else:
        context_rail.add_class("hidden")          # hide panel
        self.resource_bar.remove_class("hidden")  # SHOW resource bar
```

When context panel is open: ResourceBar has `display: none`, so "Compacting..." is **invisible**. The callback still mutates internal state on a hidden widget.

### Comparison with Agent Request UX

| Aspect | Agent Request | `/compact` Command |
|--------|--------------|-------------------|
| Loading indicator | `LoadingIndicator` shown | None |
| Viewport state | `#viewport.streaming` CSS class | None |
| Input disabled | Implicitly (event handler blocks) | Implicitly (event handler blocks) |
| Chat container output | Response streams in | Nothing written |
| Progress feedback | Spinner + streaming text | "Compacting..." in ResourceBar only |
| Post-completion | Resource bar updated | Toast notification + resource bar updated |

### Comparison with `/update` Command UX

`/update` writes progress directly to the chat container:
```python
app.chat_container.write("Checking for updates...")
app.chat_container.write(f"Installing with {pkg_mgr}...")
```

`/compact` writes **nothing** to the chat container. All feedback is either the ResourceBar label or toast notifications.

## Key Patterns / Solutions Found

- **Mutual exclusion pattern**: ResourceBar ↔ context panel via CSS `hidden` class creates blind spots for any status rendered only in the ResourceBar
- **Status callback**: Simple `Callable[[bool], None]` — only signals start/end, no progress granularity
- **No viewport state class**: Agent requests use `#viewport.streaming` for visual feedback; `/compact` has no equivalent
- **No chat container output**: Unlike `/update`, compaction writes nothing to the conversation area
- **No explicit pre-flight indicator**: The command relies on the controller calling `callback(True)` — there's no immediate visual change when the user types `/compact`

## Knowledge Gaps

1. **Exact timing**: How long does a typical compaction take? (Depends on context size + LLM latency)
2. **Error UX**: If `force_compact()` raises an exception, the `finally` block clears the status but no error toast is shown to the user
3. **Concurrent safety**: What happens if a user submits a message while `/compact` is running? The event handler awaits, so it queues, but no UI indication of this
4. **Context panel compaction display**: The context panel shows tokens/model/cost/files but has zero compaction-related information

## UX Issues Identified

### P0: No feedback when context panel is open
- **Problem**: "Compacting..." renders only in ResourceBar which is `display: none` when context panel is visible
- **Impact**: User sees nothing happening during a potentially long operation
- **Fix**: Add compaction status to context panel OR write to chat container

### P1: No immediate visual response
- **Problem**: After typing `/compact`, nothing visible happens until the controller fires `callback(True)` — could be delayed by token estimation and threshold checks
- **Fix**: Show feedback immediately in `CompactCommand.execute()` before awaiting `force_compact()`

### P2: No error feedback on exception
- **Problem**: If `force_compact()` throws, the `finally` block clears status but no toast/notification is shown
- **Fix**: Add `except Exception` with `app.notify(str(e), severity="error")`

### P3: No chat container record
- **Problem**: Compaction leaves no trace in the conversation. User can't scroll back and see that compaction happened or what it did
- **Fix**: Write a brief compaction summary to `chat_container` (matching `/update` pattern)

### P4: No viewport state feedback
- **Problem**: During agent requests, the viewport gets a `.streaming` class (visual bevel inversion). `/compact` has no equivalent, so the UI looks idle during a long operation
- **Fix**: Add `#viewport.compacting` CSS class or reuse `.streaming`

### P5: Toast notifications are ephemeral
- **Problem**: Completion/skip/error toasts auto-dismiss and are never logged to the conversation
- **Impact**: User misses results if they look away

## References

- `src/tunacode/ui/commands/compact.py` — command implementation
- `src/tunacode/core/compaction/controller.py` — controller with callback mechanism
- `src/tunacode/ui/widgets/resource_bar.py` — "Compacting..." display
- `src/tunacode/ui/app.py` — toggle logic (L419-438), status bridge (L531-556)
- `src/tunacode/ui/styles/layout.tcss` — ResourceBar.hidden CSS
- `memory-bank/research/2026-02-12_issue-382_manual-compact-boundary-map.md` — prior research on manual compact boundaries
- `memory-bank/research/2026-02-10_21-29-05_compaction-system-mapping.md` — prior compaction system mapping
