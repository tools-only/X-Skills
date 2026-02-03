---
title: "Insert Before Stream Race Condition Fix – Plan"
phase: Plan
date: "2026-01-27_17-41-42"
owner: "claude-opus-4-5"
parent_research: "memory-bank/research/2026-01-27_insert-before-stream-investigation.md"
git_commit_at_plan: "5efc5423"
tags: [plan, ui, streaming, race-condition, coding]
---

## Goal

- Fix the race condition where tool panels append to the end instead of inline when tools complete after stream cancellation or end.

## Non-Goals

- No deployment/ops changes
- No observability/metrics additions
- No new dependencies

## Scope & Assumptions

**In Scope:**
- Fix `insert_before_stream()` to correctly position tool panels even after stream ends
- Track insertion point when stream ends so late-arriving tool panels insert correctly
- Handle cancel scenario where tools complete post-cancel

**Out of Scope:**
- Tool cancellation signal propagation (separate concern)
- Throttling/batching of tool panels (LOW priority per research)

**Assumptions:**
- Textual's `mount(before=widget)` works correctly when widget exists
- Tool panels arrive via `post_message()` which is async (confirmed in research)
- Single-threaded event loop; no concurrent message handling

## Deliverables

- Modified `ChatContainer` class with insertion point tracking
- Updated `insert_before_stream()` to use tracked position
- Updated `end_stream()` and `cancel_stream()` to preserve insertion context

## Readiness

- Research complete with flow analysis and root cause identified
- Code locations verified in `chat.py` (lines 137-268)
- No external dependencies required

## Milestones

- **M1:** Insertion point tracking infrastructure
- **M2:** Fix `insert_before_stream()` logic
- **M3:** Basic test scenario (manual verification)

## Work Breakdown (Tasks)

### Task 1: Add insertion anchor tracking to ChatContainer
**Summary:** Track the last widget before which tool panels should insert, persisting after stream ends.

**Files:** `src/tunacode/ui/widgets/chat.py`

**Changes:**
1. Add `_insertion_anchor: Widget | None` attribute in `__init__`
2. In `start_stream()`: set `_insertion_anchor = None` (panels go before stream widget)
3. In `end_stream()`: capture the finalized `_current_stream` widget as `_insertion_anchor` BEFORE setting `_current_stream = None`
4. In `cancel_stream()`: capture widget reference before removal, set `_insertion_anchor` to next sibling if exists

**Acceptance Test:** After calling `end_stream()`, `_insertion_anchor` holds reference to the finalized message widget.

**Dependencies:** None
**Milestone:** M1

---

### Task 2: Update insert_before_stream() to use insertion anchor
**Summary:** Modify logic to insert before anchor when stream is None but anchor exists.

**Files:** `src/tunacode/ui/widgets/chat.py`

**Changes:**
```python
def insert_before_stream(self, renderable: RenderableType) -> None:
    widget = Static(renderable)
    widget.add_class("chat-message")

    if self._current_stream is not None:
        # Active stream: insert before streaming widget
        self.mount(widget, before=self._current_stream)
    elif self._insertion_anchor is not None:
        # Stream ended: insert before the finalized message
        self.mount(widget, before=self._insertion_anchor)
    else:
        # No context: append
        self.mount(widget)

    if self._auto_scroll:
        self.scroll_end(animate=False)
```

**Acceptance Test:** Tool panel arriving after `end_stream()` appears before the finalized agent response, not at bottom.

**Dependencies:** Task 1
**Milestone:** M2

---

### Task 3: Clear insertion anchor on new request
**Summary:** Reset anchor when starting a new stream to prevent stale positioning.

**Files:** `src/tunacode/ui/widgets/chat.py`

**Changes:**
1. In `start_stream()`: set `_insertion_anchor = None` at the start
2. In `clear()`: set `_insertion_anchor = None`

**Acceptance Test:** New request does not use stale anchor from previous request.

**Dependencies:** Task 1
**Milestone:** M2

---

### Task 4: Handle cancel scenario insertion anchor
**Summary:** On cancel, preserve insertion context for late-arriving tool panels.

**Files:** `src/tunacode/ui/widgets/chat.py`

**Changes:**
1. In `cancel_stream()`: before removing `_current_stream`, find its next sibling
2. Set `_insertion_anchor` to next sibling (or None if last child)
3. Alternative: Keep `_current_stream` hidden instead of removing? (research: simpler to just track position)

**Design Decision:** Track position via index, not sibling reference. Siblings can be fragile.

Better approach:
```python
def cancel_stream(self) -> None:
    if self._current_stream is not None:
        # Find the index of the stream widget
        children = list(self.children)
        try:
            idx = children.index(self._current_stream)
            # Anchor is widget that will be at this position after removal
            if idx > 0:
                self._insertion_anchor = children[idx - 1]
        except ValueError:
            pass
        self._current_stream.remove()
        self._current_stream = None
    self.remove_class("streaming")
```

Wait - this anchors to widget BEFORE, but we want to insert AFTER it. Let me reconsider.

**Revised approach:** Use `after=` parameter for mount instead:

```python
def cancel_stream(self) -> None:
    if self._current_stream is not None:
        children = list(self.children)
        try:
            idx = children.index(self._current_stream)
            if idx > 0:
                self._insertion_after_anchor = children[idx - 1]
        except ValueError:
            pass
        self._current_stream.remove()
        self._current_stream = None
    self.remove_class("streaming")
```

Then in `insert_before_stream()`:
```python
elif self._insertion_after_anchor is not None:
    self.mount(widget, after=self._insertion_after_anchor)
```

**Acceptance Test:** Tool panel arriving after `cancel_stream()` appears at the position where streaming was occurring.

**Dependencies:** Task 2
**Milestone:** M2

---

### Task 5: Manual verification
**Summary:** Verify fix with real tool calls in TUI.

**Files:** None (manual testing)

**Steps:**
1. Start tunacode TUI
2. Run request that triggers tool (e.g., "read file X")
3. Wait for tool panel to appear inline during streaming
4. Cancel mid-stream (Ctrl+C or escape)
5. Observe: any late tool panels should NOT appear at bottom

**Acceptance Test:** Visual confirmation that tool panels maintain correct ordering in all scenarios.

**Dependencies:** Tasks 1-4
**Milestone:** M3

## Risks & Mitigations

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Textual `after=` behavior differs from expected | Low | Medium | Verify Textual docs; test in isolation |
| Removed widget reference causes issues | Low | Low | Clear anchor when widget removed by other means |
| Stale anchor persists across sessions | Low | Low | Clear in `clear()` and `start_stream()` |

## Test Strategy

- **Task 2 acceptance:** Create unit test that mounts ChatContainer, starts stream, ends stream, then calls `insert_before_stream()` and verifies widget order.

Single test file: `tests/ui/test_chat_container.py`

## References

- Research doc: `memory-bank/research/2026-01-27_insert-before-stream-investigation.md`
- Code: `src/tunacode/ui/widgets/chat.py:244-261` (current implementation)
- Textual mount docs: https://textual.textualize.io/api/widget/#textual.widget.Widget.mount

## Final Gate

- **Plan path:** `memory-bank/plan/2026-01-27_17-41-42_insert-before-stream-fix.md`
- **Milestone count:** 3
- **Tasks ready for coding:** 5

**Next command:** `/context-engineer:execute "memory-bank/plan/2026-01-27_17-41-42_insert-before-stream-fix.md"`
