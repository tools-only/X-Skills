---
title: "Issue #311 Core Types Layer Fix – Plan"
phase: Plan
date: "2026-01-27T00:19:58"
owner: "claude-agent"
parent_research: "memory-bank/research/2026-01-27_00-10-38_issue-311-core-types-layer-violation.md"
git_commit_at_plan: "72bf42bc"
tags: [plan, architecture, types, clean-architecture]
---

## Goal

Move core-specific types from shared `types/` layer into `core/types/` to follow Clean Architecture principles where domain types live in the domain layer.

**Non-goals:**
- Changing the shared types that are correctly utils-level
- Adding re-export indirection layers
- Modifying any business logic

## Scope & Assumptions

**In scope:**
- Create `src/tunacode/core/types/` directory structure
- Move 9 core-specific types to `core/types/`
- Update 17 files in core to import from new location
- Delete 2 dead types (CommandContext, ProcessRequestCallback)
- Update DEPENDENCY_MAP.md

**Out of scope:**
- Changes to UI layer
- Changes to tools layer
- New features or behavior changes

**Assumptions:**
- The existing `core/types.py` re-export pattern can be expanded
- Tests pass after import path changes
- No circular import issues will arise

## Deliverables (DoD)

1. `src/tunacode/core/types/` directory with moved types
2. All 17 core files updated to import from `core/types/`
3. Dead code deleted from `types/dataclasses.py` and `types/callbacks.py`
4. `ruff check` passes
5. `uv run pytest` passes
6. DEPENDENCY_MAP.md updated

## Readiness (DoR)

- [x] Research document complete
- [x] Best practice verified (Clean Architecture standard)
- [x] File inventory complete (17 files, 26 imports)
- [x] Types categorized (9 to move, 2 to delete)

## Milestones

- **M1:** Create core/types/ structure and move types
- **M2:** Update all import statements
- **M3:** Delete dead code and verify

## Work Breakdown (Tasks)

### Task 1: Create core/types/ directory structure
**Summary:** Create the directory and initial `__init__.py`
**Files:** `src/tunacode/core/types/__init__.py`
**Acceptance:**
- Directory exists
- `__init__.py` exports nothing yet (placeholder)

### Task 2: Move core-specific types to core/types/
**Summary:** Move 9 types from `types/` to `core/types/`
**Types to move:**
- `ResponseState` (from dataclasses.py)
- `AgentState` (from dataclasses.py)
- `StateManagerProtocol` (from state.py)
- `SessionStateProtocol` (from state.py)
- `ConversationState` (from state_structures.py)
- `TaskState` (from state_structures.py)
- `RuntimeState` (from state_structures.py)
- `UsageState` (from state_structures.py)
- `ToolCallRegistry` (from tool_registry.py)

**Files touched:**
- Create: `src/tunacode/core/types/state.py`
- Create: `src/tunacode/core/types/protocols.py`
- Create: `src/tunacode/core/types/registry.py`
- Update: `src/tunacode/core/types/__init__.py`

**Acceptance:**
- All 9 types defined in core/types/
- Exports work via `from tunacode.core.types import X`

### Task 3: Update 17 core files to import from core/types/
**Summary:** Change import statements from `tunacode.types` to `tunacode.core.types` for moved types
**Files to update (17):**
1. `core/state.py`
2. `core/user_configuration.py`
3. `core/lsp_status.py`
4. `core/types.py` (the re-export file)
5. `core/configuration.py`
6. `core/logging/manager.py`
7. `core/agents/main.py`
8. `core/agents/agent_components/agent_config.py`
9. `core/agents/agent_components/agent_helpers.py`
10. `core/agents/agent_components/response_state.py`
11. `core/agents/agent_components/state_transition.py`
12. `core/agents/agent_components/streaming.py`
13. `core/agents/agent_components/orchestrator/orchestrator.py`
14. `core/agents/agent_components/orchestrator/tool_dispatcher.py`
15. `core/agents/agent_components/orchestrator/usage_tracker.py`
16. `core/agents/resume/sanitize.py`
17. `core/agents/resume/sanitize_debug.py`

**Acceptance:**
- All imports updated
- `ruff check` passes
- No import errors at runtime

### Task 4: Delete dead code
**Summary:** Remove unused types from shared types layer
**Delete:**
- `CommandContext` from `types/dataclasses.py`
- `ProcessRequestCallback` from `types/callbacks.py`

**Acceptance:**
- Dead code removed
- No references remain
- `ruff check` passes

### Task 5: Verify and update documentation
**Summary:** Run tests, update DEPENDENCY_MAP.md
**Files:**
- `docs/architecture/DEPENDENCY_MAP.md`

**Acceptance:**
- `uv run pytest` passes
- DEPENDENCY_MAP.md reflects new structure
- `core → types` import count reduced

## Risks & Mitigations

| Risk | Impact | Likelihood | Mitigation | Trigger |
|------|--------|------------|------------|---------|
| Circular imports | High | Medium | Careful ordering in __init__.py | Import error on startup |
| Missed import | Medium | Low | grep for old import paths | Test failure |
| Type checker errors | Medium | Medium | Run mypy, use `git commit -n` if blocked | mypy failure |

## Test Strategy

- Run existing test suite: `uv run pytest`
- Manual verification: `python -c "from tunacode.core.types import AgentState"`
- No new tests needed (structural change only)

## Tickets Created

| Ticket ID | Title | Priority | Status | Depends On |
|-----------|-------|----------|--------|------------|
| tun-3f92 | Create core/types/ directory structure | P1 | open | - |
| tun-979a | Move 9 core-specific types to core/types/ | P1 | open | tun-3f92 |
| tun-48fe | Update 17 core files to import from core/types/ | P1 | open | tun-979a |
| tun-f096 | Delete dead code (CommandContext, ProcessRequestCallback) | P2 | open | tun-48fe |
| tun-3cde | Verify tests and update DEPENDENCY_MAP.md | P2 | open | tun-48fe |

## References

- Research: `memory-bank/research/2026-01-27_00-10-38_issue-311-core-types-layer-violation.md`
- Issue: https://github.com/alchemiststudiosDOTai/tunacode/issues/311
- Clean Architecture best practice: domain types live in domain layer
