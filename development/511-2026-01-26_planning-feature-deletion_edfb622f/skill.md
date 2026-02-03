# Research – Planning Feature Deletion Map

**Date:** 2026-01-26
**Owner:** Agent
**Phase:** Research

## Goal

Map ALL layers of the planning feature in tunacode to enable clean deletion. The feature has not been actively worked on and needs to be removed.

## Findings

The planning feature is a **read-only mode** that restricts the agent to information gathering before allowing code modifications. It spans 4 layers:

### Layer 1: TOOLS (Delete Files)

| File | Purpose |
|------|---------|
| `src/tunacode/tools/present_plan.py` | Main tool implementation (46-97 lines of core logic) |
| `src/tunacode/tools/prompts/present_plan_prompt.xml` | XML prompt template for tool usage instructions |

### Layer 2: UI (Delete Files)

| File | Purpose |
|------|---------|
| `src/tunacode/ui/plan_approval.py` | Plan approval panel rendering, keyboard handling (1-3 key approval) |

### Layer 3: UI (Modify - Remove References)

| File | What to Remove |
|------|----------------|
| `src/tunacode/ui/commands/__init__.py` | `PlanCommand` class (lines ~291-321), toggles plan mode via `/plan` |
| `src/tunacode/ui/app.py` | `request_plan_approval()` method (~line 322), `_handle_plan_approval_key()` (~line 571), `pending_plan_approval` attribute (~line 101), key handling in `on_key()` (~line 546) |
| `src/tunacode/ui/repl_support.py` | `PendingPlanApprovalState` class (~line 112) |
| `src/tunacode/ui/welcome.py` | `/plan` command mention (~lines 54-55) |
| `src/tunacode/ui/widgets/status_bar.py` | `set_mode()` method - **KEEP** (shared with other modes), just verify no plan-specific logic |

### Layer 4: CORE (Modify - Remove References)

| File | What to Remove |
|------|----------------|
| `src/tunacode/core/state.py` | `plan_mode: bool` field (~line 54), `plan_approval_callback` field (~line 57) |
| `src/tunacode/core/agents/agent_components/agent_config.py` | `present_plan` tool creation and registration (~lines 433-435) |

### Layer 5: AUTHORIZATION (Modify - Simplify)

| File | What to Remove |
|------|----------------|
| `src/tunacode/tools/authorization/rules.py` | `PlanModeBlockRule` class (~lines 72-86), `PLAN_MODE_BLOCKED_TOOLS` constant (~line 12) |
| `src/tunacode/tools/authorization/context.py` | `plan_mode` field in `AuthContext` (~line 17), `from_state()` plan_mode extraction (~line 30) |
| `src/tunacode/tools/authorization/factory.py` | `PlanModeBlockRule` import and registration (~lines 6, 16) |

### Layer 6: TYPES (Modify - Remove Protocols)

| File | What to Remove |
|------|----------------|
| `src/tunacode/types/__init__.py` | Exports for planning types |
| `src/tunacode/types/state.py` | `PlanSessionProtocol` (~lines 22-27), `PlanApprovalProtocol` (~lines 29-36) |
| `src/tunacode/types/callbacks.py` | `PlanApprovalCallback` type alias (~line 64), contract docs (~lines 10-12) |

### Layer 7: CONSTANTS (Modify)

| File | What to Remove |
|------|----------------|
| `src/tunacode/constants.py` | `EXIT_PLAN_MODE_SENTINEL` (~line 105), `ToolName.PRESENT_PLAN` from `READ_ONLY_TOOLS` (~line 91) |

### Layer 8: TESTS (Delete)

| File | Test Count | Purpose |
|------|------------|---------|
| `tests/unit/core/test_present_plan.py` | 7 tests | Complete test suite for present_plan tool |

**Test functions to delete:**
- `TestPresentPlanTool.test_rejects_when_not_in_plan_mode`
- `TestPresentPlanTool.test_auto_approves_without_callback`
- `TestPresentPlanTool.test_handles_approval_callback_approve`
- `TestPresentPlanTool.test_handles_approval_callback_deny`
- `TestPresentPlanTool.test_handles_exit_sentinel`
- `TestPresentPlanTool.test_signature_preserved`
- `TestPresentPlanRegistration.test_present_plan_registered_in_agent`

### Layer 9: DOCUMENTATION (Update)

| File | Action |
|------|--------|
| `docs/codebase-map/modules/types.md` | Remove `PlanApprovalProtocol` (~line 39), `PlanApprovalCallback` (~line 60) |
| `docs/codebase-map/modules/ui-overview.md` | Remove `PlanCommand` (~line 133) |
| `docs/codebase-map/structure/02-core-directory.md` | Remove planning module mentions (~lines 33, 69) |
| `docs/codebase-map/structure/tree-structure.txt` | Remove planning module from tree (~line 67) |

### Related Tickets (Reference Only)

| File | Mention |
|------|---------|
| `.tickets/t-6670.md` | Mentions `PlanApprovalProtocol` |
| `.tickets/t-46de.md` | Mentions plan approval callbacks |
| `.tickets/t-69f0.md` | Discusses plan approval callback types |

## Recommended Deletion Order

**Phase 1: Delete standalone files (no dependents)**
1. `tests/unit/core/test_present_plan.py` - test file
2. `src/tunacode/tools/prompts/present_plan_prompt.xml` - prompt template
3. `src/tunacode/ui/plan_approval.py` - UI panel

**Phase 2: Remove from agent registration**
4. `src/tunacode/core/agents/agent_components/agent_config.py` - remove tool registration

**Phase 3: Delete tool implementation**
5. `src/tunacode/tools/present_plan.py` - tool itself

**Phase 4: Remove authorization layer references**
6. `src/tunacode/tools/authorization/factory.py` - remove rule registration
7. `src/tunacode/tools/authorization/rules.py` - remove PlanModeBlockRule class
8. `src/tunacode/tools/authorization/context.py` - remove plan_mode field

**Phase 5: Remove UI command and state**
9. `src/tunacode/ui/commands/__init__.py` - remove PlanCommand
10. `src/tunacode/ui/app.py` - remove approval handling
11. `src/tunacode/ui/repl_support.py` - remove PendingPlanApprovalState
12. `src/tunacode/ui/welcome.py` - remove /plan mention

**Phase 6: Remove core state fields**
13. `src/tunacode/core/state.py` - remove plan_mode and callback fields

**Phase 7: Remove type definitions**
14. `src/tunacode/types/callbacks.py` - remove PlanApprovalCallback
15. `src/tunacode/types/state.py` - remove plan protocols
16. `src/tunacode/types/__init__.py` - remove exports

**Phase 8: Clean constants**
17. `src/tunacode/constants.py` - remove EXIT_PLAN_MODE_SENTINEL, update READ_ONLY_TOOLS

**Phase 9: Update documentation**
18. Update all docs/codebase-map files listed above

## Known Failing Test Note

The user mentioned a failing test related to "layers" that should be **ignored** (not touched). This is likely related to the ongoing Gate 2 (dependency direction) work mentioned in CLAUDE.md. The test_present_plan.py tests should be deleted as part of this cleanup - they test planning functionality specifically.

## Key Patterns / Solutions Found

| Pattern | Description |
|---------|-------------|
| Callback Injection | UI provides `plan_approval_callback` closure to session state, tool invokes without knowing UI |
| Sentinel Value | `EXIT_PLAN_MODE_SENTINEL` distinguishes "revise plan" from "exit entirely" |
| Priority-Based Rules | `PlanModeBlockRule` at priority 100 blocks writes before other rules evaluate |
| Protocol Decoupling | `PlanApprovalProtocol` avoids circular imports between core/tools/ui |

## Knowledge Gaps

- Which specific test is "failing due to layers" that should be ignored?
- Any external integrations expecting `/plan` command?

## References

- `src/tunacode/tools/present_plan.py` - Main tool implementation
- `src/tunacode/ui/plan_approval.py` - UI panel
- `src/tunacode/tools/authorization/rules.py` - Authorization blocking
- `tests/unit/core/test_present_plan.py` - All planning tests
