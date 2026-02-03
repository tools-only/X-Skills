---
title: "Planning Feature Deletion – Plan"
phase: Plan
date: "2026-01-26T14:30:00Z"
owner: "agent"
parent_research: "memory-bank/research/2026-01-26_planning-feature-deletion.md"
git_commit_at_plan: "ba2ac938"
tags: [plan, deletion, cleanup, planning-feature]
---

## Goal

**Remove the planning feature entirely from tunacode.** The feature (a read-only mode restricting agent to information gathering before code modifications) is not actively used and adds maintenance burden across 9 layers.

**Non-goals:**
- Refactoring related code beyond deletion
- Adding new features to replace planning
- Touching unrelated failing tests (Gate 2 dependency work)

## Scope & Assumptions

**In scope:**
- Delete 4 files (tool, prompt, UI panel, tests)
- Modify 17 files to remove planning references
- Update 4 documentation files

**Out of scope:**
- The failing "layers" test (Gate 2 work, unrelated)
- Any functionality enhancements

**Assumptions:**
- No external integrations depend on `/plan` command
- All planning code is captured in the research doc
- Branch `delete-plan` is the working branch

## Deliverables (DoD)

| Artifact | Acceptance Criteria |
|----------|---------------------|
| Clean deletion | All 4 files deleted, no orphaned imports |
| No regressions | `uv run pytest` passes (except known Gate 2 failures) |
| No type errors | `uv run mypy` shows no NEW errors |
| Clean lint | `uv run ruff check .` passes |
| Updated docs | All codebase-map references removed |

## Readiness (DoR)

- [x] Research document complete with full file/line mapping
- [x] Git branch `delete-plan` exists
- [x] Working directory clean

## Milestones

- **M1:** Delete standalone files (tests, prompts, UI panel)
- **M2:** Remove tool registration and delete tool implementation
- **M3:** Remove authorization layer references
- **M4:** Remove UI command and state references
- **M5:** Remove types, constants, and update documentation

## Work Breakdown (Tasks)

### Task 1: Delete Standalone Files (M1)
**Summary:** Delete files with no dependents
**Owner:** agent
**Dependencies:** none
**Files:**
- DELETE `tests/unit/core/test_present_plan.py`
- DELETE `src/tunacode/tools/prompts/present_plan_prompt.xml`
- DELETE `src/tunacode/ui/plan_approval.py`

**Acceptance Tests:**
- Files no longer exist
- No import errors when running pytest

---

### Task 2: Remove Tool Registration & Delete Tool (M2)
**Summary:** Unregister from agent, then delete tool
**Owner:** agent
**Dependencies:** Task 1
**Files:**
- MODIFY `src/tunacode/core/agents/agent_components/agent_config.py` - remove present_plan registration (~lines 433-435)
- DELETE `src/tunacode/tools/present_plan.py`

**Acceptance Tests:**
- Agent initializes without present_plan tool
- No import errors

---

### Task 3: Remove Authorization Layer (M3)
**Summary:** Remove plan mode blocking rules
**Owner:** agent
**Dependencies:** Task 2
**Files:**
- MODIFY `src/tunacode/tools/authorization/factory.py` - remove PlanModeBlockRule import/registration
- MODIFY `src/tunacode/tools/authorization/rules.py` - remove PlanModeBlockRule class, PLAN_MODE_BLOCKED_TOOLS
- MODIFY `src/tunacode/tools/authorization/context.py` - remove plan_mode field

**Acceptance Tests:**
- Authorization system works without plan_mode
- Tool authorization tests pass

---

### Task 4: Remove UI Command & State (M4)
**Summary:** Remove /plan command and approval handling
**Owner:** agent
**Dependencies:** Task 3
**Files:**
- MODIFY `src/tunacode/ui/commands/__init__.py` - remove PlanCommand class
- MODIFY `src/tunacode/ui/app.py` - remove request_plan_approval(), _handle_plan_approval_key(), pending_plan_approval
- MODIFY `src/tunacode/ui/repl_support.py` - remove PendingPlanApprovalState
- MODIFY `src/tunacode/ui/welcome.py` - remove /plan mention
- MODIFY `src/tunacode/core/state.py` - remove plan_mode, plan_approval_callback

**Acceptance Tests:**
- App starts without plan-related attributes
- /plan command no longer available
- Welcome message doesn't mention /plan

---

### Task 5: Remove Types, Constants & Update Docs (M5)
**Summary:** Clean up type definitions and documentation
**Owner:** agent
**Dependencies:** Task 4
**Files:**
- MODIFY `src/tunacode/types/callbacks.py` - remove PlanApprovalCallback
- MODIFY `src/tunacode/types/state.py` - remove PlanSessionProtocol, PlanApprovalProtocol
- MODIFY `src/tunacode/types/__init__.py` - remove planning exports
- MODIFY `src/tunacode/constants.py` - remove EXIT_PLAN_MODE_SENTINEL, update READ_ONLY_TOOLS
- MODIFY `docs/codebase-map/modules/types.md`
- MODIFY `docs/codebase-map/modules/ui-overview.md`
- MODIFY `docs/codebase-map/structure/02-core-directory.md`
- MODIFY `docs/codebase-map/structure/tree-structure.txt`

**Acceptance Tests:**
- No orphaned type exports
- Documentation accurate
- Full test suite passes

---

## Risks & Mitigations

| Risk | Impact | Likelihood | Mitigation | Trigger |
|------|--------|------------|------------|---------|
| Hidden dependency on planning | High | Low | Research doc is comprehensive | Import error during deletion |
| Merge conflicts with master | Medium | Low | Working on dedicated branch | Rebase before merge |
| Type errors cascade | Medium | Medium | Delete in dependency order | Mypy fails after deletion |

## Test Strategy

- Run `uv run pytest` after each task to catch regressions immediately
- One verification run of full test suite at end
- No new tests needed (deleting feature, not adding)

## References

- Research doc: `memory-bank/research/2026-01-26_planning-feature-deletion.md`
- Main tool: `src/tunacode/tools/present_plan.py`
- UI panel: `src/tunacode/ui/plan_approval.py`
- Authorization: `src/tunacode/tools/authorization/rules.py`

## Tickets Created (max 5)

| Ticket ID | Title | Priority | Status |
|-----------|-------|----------|--------|
| tun-2ef4 | Delete standalone planning files (tests, prompts, UI panel) | 1 | open |
| tun-6d0a | Remove present_plan tool registration and delete tool | 1 | open |
| tun-0de7 | Remove authorization layer planning references | 2 | open |
| tun-216b | Remove UI command and state references | 2 | open |
| tun-4014 | Remove planning types, constants, and update docs | 2 | open |

## Dependencies

```
tun-2ef4 (Phase 1: Delete standalone files)
    └── tun-6d0a (Phase 2: Remove tool registration)
        └── tun-0de7 (Phase 3: Authorization layer)
            └── tun-216b (Phase 4: UI command/state)
                └── tun-4014 (Phase 5: Types/constants/docs)
```
