---
title: Centralized Panel Width Handling
link: panel-width-refactoring
type: delta
path: src/tunacode/ui/renderers/panel_widths.py
depth: 0
seams: [M]
ontological_relations:
  - relates_to: [[tool-renderers]]
  - affects: [[ui-overview]]
  - fixes: [[panel-width-indirection]]
tags:
  - refactor
  - ui
  - tool-renderers
  - gate-5
created_at: 2026-01-16T16:51:33Z
updated_at: 2026-01-16T16:51:33Z
uuid: e5ac0f53-c6c5-4536-a4d7-d6d7488825bf
---

# Centralized Panel Width Handling

## Summary

PR #244 centralized panel width calculations into a single module (`panel_widths.py`) and removed redundant width logic from 9+ tool renderers. This refactoring applies Gate 5 (Indirection Requires Verification) by replacing `expand=True` indirection with explicit `width=` parameters.

## Context

- File: `src/tunacode/ui/renderers/panel_widths.py`
- PR: #244
- Commit: `f0bfbb1`
- Surfaces: Tool panels now have consistent, verifiable widths

## Root Cause

**Indirection hiding actual behavior.**

Tool renderers were using `Panel(content, expand=True)` and trusting that panels would fill the terminal width. But:

1. `expand=True` expands to `scrollable_content_region`, not terminal width
2. RichLog's internal padding/gutter subtracted 4-8 chars
3. Each renderer had its own width calculation logic
4. No single place controlled the final panel width

This violated Gate 5: we were passing a wish (`expand=True`) instead of verifying the output.

```
BEFORE (indirect - unverifiable):
┌─────────────────────────────────────────────────────────────┐
│ Panel(content, expand=True)                                 │
│   → RichLog.write(panel)                                    │
│   → RichLog decides width based on scrollable_content_region│
│   → Panel renders narrower than expected                    │
│   → No way to verify what width was actually used           │
└─────────────────────────────────────────────────────────────┘

AFTER (direct - verifiable):
┌─────────────────────────────────────────────────────────────┐
│ frame_width = tool_panel_frame_width(max_line_width)        │
│ Panel(content, width=frame_width)                           │
│   → Width is explicit, testable, traceable                  │
│   → panel_widths.py is the single source of truth           │
└─────────────────────────────────────────────────────────────┘
```

## Changes

### New Module: `panel_widths.py`

```python
def tool_panel_frame_width(max_line_width: int) -> int:
    """Return tool panel frame width including borders and padding.

    Preconditions:
        - max_line_width is already clamped to a positive value.
    """
    return max_line_width + TOOL_PANEL_HORIZONTAL_INSET
```

### Simplified Renderers

9 tool renderers were simplified by removing local width calculations:

| Renderer | Lines Removed | Change |
|----------|---------------|--------|
| `glob.py` | ~20 | Removed width calculation |
| `grep.py` | ~15 | Removed width calculation |
| `read_file.py` | ~10 | Uses passed `max_line_width` |
| `update_file.py` | ~15 | Uses passed `max_line_width` |
| `write_file.py` | ~5 | Uses passed `max_line_width` |
| `list_dir.py` | ~5 | Uses passed `max_line_width` |
| `bash.py` | ~3 | Uses passed `max_line_width` |
| `research.py` | ~10 | Uses passed `max_line_width` |
| `web_fetch.py` | ~3 | Uses passed `max_line_width` |

### Base Renderer Changes

`base.py` now consistently passes `max_line_width` through the render chain:

```python
def render(self, args, result, duration_ms, max_line_width) -> RenderableType | None:
    # max_line_width flows to build_viewport() and build_status()
```

### Truncation Logic

`panels.py` simplified `_truncate_content()` - it now trusts the width passed to it rather than computing its own.

## Behavioral Impact

- Tool panels render at consistent, predictable widths
- Width is verifiable: `max_line_width + TOOL_PANEL_HORIZONTAL_INSET`
- Single place to change width policy (`panel_widths.py`)
- Renderers are simpler, focused on content not layout math

## Gate 5 Application

This refactoring demonstrates Gate 5: **Indirection Requires Verification**.

| Before | After |
|--------|-------|
| `expand=True` - who decides? | `width=frame_width` - explicit |
| 9 different width calculations | 1 function in `panel_widths.py` |
| Can't verify actual width | Width is a traceable value |
| Trust the framework | Verify the output |

**Rule:** When delegating a decision (like "expand to fill"), verify the OUTPUT, not just the INPUT. If you can't observe the final value, you're guessing.

## Files Changed

- `src/tunacode/ui/renderers/panel_widths.py` (new)
- `src/tunacode/ui/renderers/tools/*.py` (9 renderers simplified)
- `src/tunacode/ui/renderers/panels.py` (truncation logic)
- `src/tunacode/ui/renderers/tools/base.py` (hook-arrow pattern)
- `src/tunacode/ui/app.py` (cleanup)
- `tests/*.py` (updated)

## Related Cards

- [[tool-renderers]] - Tool renderer architecture
- [[ui-overview]] - UI module documentation
- [[nextstep_panels]] - NeXTSTEP design principles
