---
title: "Issue #314 Layer 2 Lateral Coupling Fix – Plan"
phase: Plan
date: "2026-01-27T21:30:00"
owner: "claude-agent"
parent_research: "memory-bank/research/2026-01-27_13-11-55_issue-314-lateral-coupling-tools-indexing-lsp.md"
git_commit_at_plan: "c5fe6f95"
tags: [plan, architecture, layer-2, lateral-coupling]
---

## Goal

Eliminate all lateral coupling violations between Layer 2 peer modules (`tools`, `indexing`, `lsp`, `infrastructure`). These modules should only depend downward to foundation layers (`utils`, `types`, `configuration`), never sideways to each other.

**Non-goals:**
- Major architectural changes beyond fixing the 3 identified violations
- Adding new features to affected tools
- Changing tool behavior (only dependency flow)

## Scope & Assumptions

**In scope:**
- 3 lateral imports to eliminate:
  1. `tools/glob.py:12` → `indexing.CodeIndex`
  2. `tools/decorators.py:73` → `lsp.{get_diagnostics, format_diagnostics}`
  3. `tools/lsp_status.py:5` → `lsp.servers.get_server_command`

**Out of scope:**
- Other Layer 2 modules (`infrastructure`)
- UI layer changes
- New test files (existing tests must pass)

**Assumptions:**
- PR #317 has merged (core→utils violations fixed)
- No new lateral violations added since research
- LSP and indexing modules remain stable during refactor

## Deliverables (DoD)

1. **Zero lateral imports between Layer 2 modules**
   - `grep -rn "from tunacode.indexing" src/tunacode/tools/` returns empty
   - `grep -rn "from tunacode.lsp" src/tunacode/tools/` returns empty

2. **All existing tests pass**
   - `uv run pytest` exits 0

3. **Type checking passes**
   - `uv run mypy src/tunacode/` no new errors

4. **Architecture documentation updated**
   - `docs/architecture/NEW_layers.html` reflects resolved state

## Readiness (DoR)

- [x] Research document complete
- [x] Violations verified in current codebase
- [x] PR #317 merged (prerequisite)
- [x] Architecture diagram ready

## Milestones

### M1: Fix Violation 3 (lsp_status.py)
**Effort:** Low | **Risk:** Low

Move `get_server_command()` to configuration layer OR merge `tools/lsp_status.py` into `core/lsp_status.py`.

Recommended: **Merge into core** - eliminates unnecessary indirection.

### M2: Fix Violation 2 (decorators.py → lsp)
**Effort:** Medium | **Risk:** Medium

Remove LSP orchestration from `@file_tool` decorator. Core should handle post-write diagnostics, not the decorator.

Recommended: **Callback injection** - decorator accepts optional `on_write_callback`, core injects LSP check.

### M3: Fix Violation 1 (glob.py → indexing)
**Effort:** Medium-High | **Risk:** Medium

Two viable options:
- **Option A:** Protocol abstraction - define `FileRegistry` in `types/`, inject from core
- **Option C:** Remove indexing check entirely from glob (simplest)

Recommended: **Option C** - glob should not need to know about indexing state. The performance optimization can be handled at core level if needed.

### M4: Validation & Documentation
**Effort:** Low | **Risk:** None

Run validation commands, update architecture docs, close issue.

## Work Breakdown (Tasks)

| ID | Task | Milestone | Dependencies | Files |
|----|------|-----------|--------------|-------|
| T1 | Merge tools/lsp_status.py into core/lsp_status.py | M1 | None | `tools/lsp_status.py`, `core/lsp_status.py` |
| T2 | Remove LSP import from decorators.py, add callback parameter | M2 | T1 | `tools/decorators.py`, `core/agents/main.py` |
| T3 | Remove CodeIndex import from glob.py | M3 | T2 | `tools/glob.py` |
| T4 | Validate no lateral imports, run tests, update docs | M4 | T3 | `docs/architecture/NEW_layers.html` |

### Task Details

**T1: Merge lsp_status into core**
- Move `get_lsp_status()` logic from `tools/lsp_status.py` to `core/lsp_status.py`
- Delete `tools/lsp_status.py`
- Update any imports
- Acceptance: `tools/lsp_status.py` deleted, no `tools→lsp` import for status

**T2: Remove LSP from decorators**
- Add optional `diagnostics_callback: Callable[[Path], str | None] = None` to `@file_tool`
- Remove lazy import of `tunacode.lsp` in `_get_lsp_diagnostics()`
- Core injects the callback when registering tools
- Acceptance: `grep "from tunacode.lsp" decorators.py` returns empty

**T3: Remove CodeIndex from glob**
- Remove `from tunacode.indexing import CodeIndex` line
- Remove `_get_code_index()` function and `_glob_with_index()` optimization
- Keep only filesystem scan path
- Acceptance: `grep "from tunacode.indexing" glob.py` returns empty

**T4: Validation**
- Run lateral import checks
- Run `uv run pytest`
- Update `NEW_layers.html` to show 0 violations
- Close issue #314

## Risks & Mitigations

| Risk | Impact | Likelihood | Mitigation |
|------|--------|------------|------------|
| Glob performance degrades without index | Medium | Low | Index was optimization, not critical path |
| LSP diagnostics break | Medium | Low | Callback pattern well-understood, tests exist |
| Core becomes too large | Low | Low | Minimal code movement, just orchestration |

## Test Strategy

Existing tests must pass. No new tests required for this refactor since:
- We're simplifying, not adding features
- Tests already cover tool behavior
- Integration tests validate end-to-end

## References

- Research: `memory-bank/research/2026-01-27_13-11-55_issue-314-lateral-coupling-tools-indexing-lsp.md`
- Issue: https://github.com/alchemiststudiosDOTai/tunacode/issues/314
- Architecture: `docs/architecture/NEW_layers.html`
- Prior PR #317: commit `076cbf3c`

## Tickets Created

| Ticket ID | Title | Priority | Status |
|-----------|-------|----------|--------|
| tun-2ddf | Merge tools/lsp_status.py into core/lsp_status.py | P1 | open |
| tun-be3d | Remove LSP import from file_tool decorator | P2 | open |
| tun-5681 | Remove CodeIndex import from glob.py | P2 | open |
| tun-1012 | Validate lateral coupling fix and close issue #314 | P3 | open |

## Dependencies

```
tun-2ddf → tun-be3d → tun-5681 → tun-1012
```

Sequential execution - each task builds on the previous.
