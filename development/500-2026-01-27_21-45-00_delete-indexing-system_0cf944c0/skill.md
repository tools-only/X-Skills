---
title: "Delete Indexing System – Plan"
phase: Plan
date: "2026-01-27T21:45:00"
owner: "claude-agent"
parent_research: "memory-bank/research/2026-01-27_13-11-55_issue-314-lateral-coupling-tools-indexing-lsp.md"
git_commit_at_plan: "c5fe6f95"
tags: [plan, architecture, indexing, deletion]
---

## Goal

Delete the entire indexing system. It's premature optimization that adds startup latency for marginal glob improvement.

**Rationale:**
- Only ONE tool (glob.py) uses the index
- Every startup: scan filesystem to count → scan again to build → hold paths in memory
- The "benefit": iterate a cached list instead of os.scandir()
- This is startup tax for negligible gain

**Non-goals:**
- Re-implementing indexing elsewhere
- Adding new optimization to glob
- Changing glob's functional behavior

## Scope & Assumptions

**DELETE (Production Code):**

| File/Dir | Lines | Purpose |
|----------|-------|---------|
| `src/tunacode/indexing/` | ~490 | Entire module |
| `src/tunacode/core/indexing_service.py` | 137 | Orchestration wrapper |
| `docs/codebase-map/modules/indexing.md` | 97 | Module docs |
| **Total** | **~724** | |

**MODIFY (Production Code):**

| File | Changes |
|------|---------|
| `tools/glob.py` | Remove CodeIndex import, delete `_get_code_index()`, delete `_glob_with_index()` |
| `core/state.py` | Remove `indexing_service` property |
| `ui/app.py` | Delete `_run_startup_index()` method and its call |

**MODIFY (Tests):**

| File | Changes |
|------|---------|
| `tests/test_dependency_layers.py` | Remove "indexing" from LAYERS, allowed imports |
| `tests/architecture/test_import_order.py` | Remove "indexing" from SHARED_LAYER_MODULES |

**MODIFY (Documentation):**

| File | Changes |
|------|---------|
| `docs/codebase-map/MAP.md` | Remove indexing entries (3 lines) |
| `docs/architecture/DEPENDENCY_MAP.md` | Remove indexing layer and rows |
| `docs/architecture/NEW_layers.html` | Remove indexing violation references |
| `docs/codebase-map/modules/tools-overview.md` | Remove CodeIndex line |

**Assumptions:**
- glob.py already has filesystem fallback (confirmed)
- No other code paths depend on CodeIndex (confirmed by research)

## Deliverables (DoD)

1. **Zero indexing imports anywhere**
   - `grep -r "from tunacode.indexing" src/tunacode/` returns empty
   - `grep -r "CodeIndex" src/tunacode/` returns empty
   - `grep -r "IndexingService" src/tunacode/` returns empty

2. **No startup indexing**
   - App starts without "Code cache built" message
   - No filesystem scans at startup

3. **Glob still works**
   - `uv run pytest` passes
   - Manual glob test succeeds

4. **Clean architecture**
   - tools → indexing violation eliminated (issue #314)

## Milestones

### M1: Delete indexing/ module and core wrapper
**Effort:** Low | **Files:** 4

1. Delete `src/tunacode/indexing/` directory
2. Delete `src/tunacode/core/indexing_service.py`

### M2: Update consumers (glob, state, ui)
**Effort:** Medium | **Files:** 3

1. Clean `tools/glob.py` - remove index optimization
2. Clean `core/state.py` - remove indexing_service property
3. Clean `ui/app.py` - remove startup indexing

### M3: Update tests and docs
**Effort:** Low | **Files:** 6+

1. Update layer tests
2. Delete indexing.md
3. Update DEPENDENCY_MAP.md, MAP.md, NEW_layers.html, tools-overview.md

### M4: Validation
**Effort:** Low

Run all verification commands, close issue #314.

## Work Breakdown (Tasks)

| ID | Task | Milestone | Files |
|----|------|-----------|-------|
| T1 | Delete indexing/ dir and core/indexing_service.py | M1 | 4 files |
| T2 | Clean glob.py, state.py, ui/app.py | M2 | 3 files |
| T3 | Update tests and docs | M3 | 6+ files |
| T4 | Validate and close issue #314 | M4 | - |

### Task Details

**T1: Delete indexing module**
- `rm -rf src/tunacode/indexing/`
- `rm src/tunacode/core/indexing_service.py`
- Acceptance: directories/files gone

**T2: Clean consumers**
- `tools/glob.py`:
  - Remove line 12: `from tunacode.indexing import CodeIndex`
  - Remove lines 67-77: index lookup logic
  - Delete `_get_code_index()` function (lines 126-135)
  - Delete `_glob_with_index()` function (lines 176-200)

- `core/state.py`:
  - Remove TYPE_CHECKING import of IndexingService
  - Remove `_indexing_service` attribute
  - Delete `indexing_service` property

- `ui/app.py`:
  - Delete `_run_startup_index()` method
  - Remove call to it from startup sequence

- Acceptance: No import errors, glob uses filesystem only

**T3: Update tests and docs**
- `tests/test_dependency_layers.py`: Remove indexing from LAYERS
- `tests/architecture/test_import_order.py`: Remove indexing from SHARED_LAYER_MODULES
- Delete `docs/codebase-map/modules/indexing.md`
- Update MAP.md, DEPENDENCY_MAP.md, NEW_layers.html, tools-overview.md

- Acceptance: Tests pass, no indexing references in docs

**T4: Validate**
```bash
# No indexing imports
grep -r "from tunacode.indexing" src/tunacode/
grep -r "CodeIndex" src/tunacode/
grep -r "IndexingService" src/tunacode/

# Tests pass
uv run pytest

# No doc references
grep -r "indexing" docs/ | grep -v ".git"
```

## Risks & Mitigations

| Risk | Impact | Likelihood | Mitigation |
|------|--------|------------|------------|
| Glob slower on large repos | Low | Medium | os.scandir is fast, acceptable trade-off |
| Missed reference | Low | Low | Research was thorough, run grep after |

## Test Strategy

Existing tests must pass. No new tests - we're deleting code, not adding features.

## References

- Issue: https://github.com/alchemiststudiosDOTai/tunacode/issues/314
- Research: `memory-bank/research/2026-01-27_13-11-55_issue-314-lateral-coupling-tools-indexing-lsp.md`
- Prior plan (superseded): `memory-bank/plan/2026-01-27_21-30-00_issue-314-lateral-coupling.md`

## Tickets

| Ticket ID | Title | Priority | Status |
|-----------|-------|----------|--------|
| tun-27ea | Delete indexing/ module and core/indexing_service.py | P1 | open |
| tun-5f58 | Clean glob.py, state.py, ui/app.py | P1 | open |
| tun-8218 | Update tests and docs for indexing removal | P2 | open |

## Dependencies

```
tun-27ea → tun-5f58 → tun-8218 → tun-2ddf (LSP fix)
```

## Notes on Issue #314

This plan fully resolves issue #314 Violation 1 (tools/glob.py → indexing) by deleting the entire indexing system rather than refactoring around it.

The other two violations (Violations 2 and 3) remain separate tickets:
- tun-2ddf: Merge tools/lsp_status.py into core
- tun-be3d: Remove LSP import from decorators
