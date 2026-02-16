# Requirements: Finish troubleshoot→debug Rename

**Feature**: finish-rename
**Status**: APPROVED
**Created**: 2026-01-30

## Context

The z-shortcut-and-rename feature was partially completed. All `zerg/` package code (commands, diagnostics, CLI, tests) was renamed. Three areas remain:

1. `.zerg/` runtime scripts (troubleshoot.py, test_troubleshoot.py)
2. Installed commands in `.claude/commands/` (stale `zerg:troubleshoot.md`)
3. Documentation/backlog references to old paths

## Functional Requirements

### FR-1: Rename .zerg/ runtime scripts
- `.zerg/troubleshoot.py` → `.zerg/debug.py` with all class/docstring renames
- `.zerg/tests/test_troubleshoot.py` → `.zerg/tests/test_debug.py` with import updates

### FR-2: Fix installed commands
- Reinstall commands so `.claude/commands/zerg:debug.md` replaces `zerg:troubleshoot.md`
- Verify z: shortcuts also updated

### FR-3: Update documentation references
- `tasks/COVERAGE-100-BACKLOG.md` — update old file path references
- `claudedocs/backlog.md` — mark item 11 complete
- `.gsd/tasks/` historical files — update command/file references
- `zerg/data/commands/zerg:debug.md` — fix one stray "troubleshooting" ref

### FR-4: Quality gates
- No stray "troubleshoot" references outside historical/English-usage contexts
- ruff + mypy clean
- Tests pass

## Non-Goals
- Changing README.md/PROJECT_INSTRUCTIONS.md generic English uses of "troubleshoot" (e.g., "good for troubleshooting")
- Modifying `.gsd/specs/troubleshoot-enhancement/` (historical spec, preserved as-is)
- Modifying `.gsd/specs/z-shortcut-and-rename/` (historical spec, preserved as-is)
