---
title: Resolve four mypy contract violations with source-level typing fixes
type: delta
link: mypy-type-contract-fixes
path: src/tunacode
depth: 1
seams: [S, M]
ontological_relations:
  - affects: [[ui]]
  - affects: [[core]]
  - affects: [[tools]]
  - affects: [[typing-contracts]]
tags:
  - mypy
  - typing
  - literals
  - rich
  - error-handling
created_at: 2026-02-22T15:12:48-06:00
updated_at: 2026-02-22T15:12:48-06:00
uuid: febc9ce8-4061-4ee8-b2a1-770f0997a7ca
---

# Resolve four mypy contract violations with source-level typing fixes

## Summary

Fixed the reported mypy errors by tightening type contracts at their source:

1. `editor.py`: added explicit Rich style narrowing so `Text.stylize()` always receives `StyleType`.
2. `web_fetch.py`: marked `_handle_http_error()` as `NoReturn` because it always raises.
3. `compaction/types.py`: made compaction status constants literal-typed (`Final[Literal[...]]`) and defined `CompactionStatus` as a type alias.
4. `shell_runner.py`: unpacked `render_bash()` tool render tuple and returned only the renderable content as required by `RenderableType`.

No shims, no fallback branches, and no cast-based suppression were introduced.

## Files changed

- `src/tunacode/ui/widgets/editor.py`
- `src/tunacode/tools/web_fetch.py`
- `src/tunacode/core/compaction/types.py`
- `src/tunacode/ui/shell_runner.py`

## Verification

- `uv run mypy src/tunacode/ui/widgets/editor.py src/tunacode/tools/web_fetch.py src/tunacode/core/compaction/controller.py src/tunacode/ui/shell_runner.py src/tunacode/core/compaction/types.py` ✅
- `uv run ruff check src/tunacode/ui/widgets/editor.py src/tunacode/tools/web_fetch.py src/tunacode/core/compaction/types.py src/tunacode/ui/shell_runner.py` ✅
