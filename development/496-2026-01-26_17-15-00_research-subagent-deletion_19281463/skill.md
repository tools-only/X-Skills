---
title: "Research Sub-Agent Feature Deletion – Plan"
phase: Plan
date: "2026-01-26T17:15:00"
owner: "agent"
parent_research: "memory-bank/research/2026-01-26_16-41-59_research-subagent-deletion.md"
git_commit_at_plan: "b6eb3990"
tags: [plan, deletion, cleanup, research-agent]
---

## Goal

Delete the disabled `research_codebase` sub-agent feature completely from the codebase. This is a cleanup task—the feature is already disabled in production (commented out in agent_config.py).

**Non-goals:**
- Re-enabling or fixing the research feature
- Adding replacement functionality
- Refactoring adjacent code

## Scope & Assumptions

**In scope:**
- Delete 7 files (2 core, 1 UI renderer, 4 prompt sections)
- Modify 17 files to remove references
- Update documentation to match reality

**Out of scope:**
- Any new features
- Test additions beyond verification

**Assumptions:**
- Feature is confirmed disabled (verified: agent_config.py:37-39 commented)
- No external dependencies on these modules
- Follows same pattern as planning feature deletion (completed successfully)

## Deliverables (DoD)

1. All 7 files deleted
2. All 17 files modified with references removed
3. `ruff check` passes
4. `pytest` passes (no import errors)
5. No remaining references to `research_codebase`, `RESEARCH_TEMPLATE`, `ResearchResult`

## Readiness (DoR)

- [x] Research document complete
- [x] Files to delete verified exist
- [x] Git state clean (only untracked research doc)
- [x] Similar deletion pattern documented (planning feature)

## Milestones

- **M1**: Delete standalone files (prompts, renderer)
- **M2**: Delete agent core files
- **M3**: Remove constants, config, template references
- **M4**: Update tests and documentation
- **M5**: Final verification (ruff, pytest, grep)

## Work Breakdown (Tasks)

### Task 1: Delete Prompt Sections Directory
**Summary:** Remove `src/tunacode/prompts/research/` directory (4 XML files)
**Owner:** agent
**Dependencies:** None
**Milestone:** M1
**Files:**
- DELETE `src/tunacode/prompts/research/sections/agent_role.xml`
- DELETE `src/tunacode/prompts/research/sections/tool_use.xml`
- DELETE `src/tunacode/prompts/research/sections/constraints.xml`
- DELETE `src/tunacode/prompts/research/sections/output_format.xml`
- DELETE `src/tunacode/prompts/research/sections/` directory
- DELETE `src/tunacode/prompts/research/` directory
**Acceptance Tests:**
- Directory `src/tunacode/prompts/research/` does not exist
- `grep -r "research" src/tunacode/prompts/` returns no results

### Task 2: Delete UI Renderer and Update Registry
**Summary:** Remove research renderer and clean up imports/registry
**Owner:** agent
**Dependencies:** None
**Milestone:** M1
**Files:**
- DELETE `src/tunacode/ui/renderers/tools/research.py`
- MODIFY `src/tunacode/ui/renderers/tools/__init__.py` (remove import, export)
- MODIFY `src/tunacode/ui/renderers/panels.py` (remove registry entry)
**Acceptance Tests:**
- File `src/tunacode/ui/renderers/tools/research.py` does not exist
- `from tunacode.ui.renderers.tools import render_research_codebase` raises ImportError
- No `research` in renderer registry

### Task 3: Delete Agent Core Files
**Summary:** Remove delegation_tools.py and research_agent.py
**Owner:** agent
**Dependencies:** Task 1, Task 2
**Milestone:** M2
**Files:**
- DELETE `src/tunacode/core/agents/delegation_tools.py`
- DELETE `src/tunacode/core/agents/research_agent.py`
- MODIFY `src/tunacode/core/agents/agent_components/agent_config.py` (remove commented import)
**Acceptance Tests:**
- Both files deleted
- `from tunacode.core.agents import delegation_tools` raises ImportError
- `from tunacode.core.agents import research_agent` raises ImportError

### Task 4: Remove Constants, Config, and Template References
**Summary:** Clean up constants, settings, dispatcher, and prompting template
**Owner:** agent
**Dependencies:** Task 3
**Milestone:** M3
**Files:**
- MODIFY `src/tunacode/constants.py` (remove ToolName.RESEARCH_CODEBASE, READ_ONLY_TOOLS entry)
- MODIFY `src/tunacode/configuration/settings.py` (remove from tool list)
- MODIFY `src/tunacode/core/agents/agent_components/orchestrator/tool_dispatcher.py` (remove constants)
- MODIFY `src/tunacode/core/agents/agent_components/agent_helpers.py` (remove research case)
- MODIFY `src/tunacode/core/prompting/templates.py` (remove RESEARCH_TEMPLATE)
- MODIFY `src/tunacode/core/prompting/__init__.py` (remove import/export)
**Acceptance Tests:**
- `ToolName.RESEARCH_CODEBASE` raises AttributeError
- `RESEARCH_TEMPLATE` not exported from prompting module
- `ruff check` passes

### Task 5: Update Tests, Scripts, and Documentation
**Summary:** Remove test cases and update docs to match reality
**Owner:** agent
**Dependencies:** Task 4
**Milestone:** M4
**Files:**
- MODIFY `tests/integration/tools/test_tool_dispatcher_coverage.py`
- MODIFY `tests/integration/core/test_tool_call_lifecycle.py`
- MODIFY `scripts/preview_tool_panels.py`
- MODIFY `docs/codebase-map/modules/core-agents.md`
- MODIFY `docs/codebase-map/modules/core-prompting.md`
- MODIFY `docs/codebase-map/modules/ui-overview.md`
- MODIFY `docs/codebase-map/structure/02-core-directory.md`
- MODIFY `docs/codebase-map/structure/tree-structure.txt`
**Acceptance Tests:**
- `pytest` passes
- `grep -r "research_codebase" docs/` returns no results
- `grep -r "RESEARCH_TEMPLATE" docs/` returns no results

## Risks & Mitigations

| Risk | Impact | Likelihood | Mitigation | Trigger |
|------|--------|------------|------------|---------|
| Hidden import dependency | Medium | Low | Grep for all imports before deletion | Import error in pytest |
| Test failure from missing mock | Medium | Low | Run pytest after each deletion phase | Test failure |
| Documentation drift | Low | Low | Update docs in same commit as code | PR review |

## Test Strategy

- **Verification only**: Run existing tests after each milestone
- No new tests needed (deleting dead code)
- Final verification: `ruff check . && pytest`

## References

- Research doc: `memory-bank/research/2026-01-26_16-41-59_research-subagent-deletion.md`
- Similar pattern: `memory-bank/plan/2026-01-26_14-30-00_planning-feature-deletion.md`
- Bug that disabled feature: `memory-bank/research/2026-01-25_limited-read-file-parameter-mismatch.md`

## Tickets Created

| Ticket ID | Title | Priority | Status |
|-----------|-------|----------|--------|
| tun-8d48 | Delete research prompt sections directory | P2 | open |
| tun-a63b | Delete research UI renderer and update registry | P2 | open |
| tun-5a52 | Delete delegation_tools.py and research_agent.py | P2 | open |
| tun-c36f | Remove research constants, config, and template refs | P2 | open |
| tun-be14 | Update tests, scripts, and documentation | P3 | open |

## Dependencies

```
tun-8d48 ──┐
           ├──> tun-5a52 ──> tun-c36f ──> tun-be14
tun-a63b ──┘
```

- `tun-5a52` depends on `tun-8d48` and `tun-a63b` (can start after prompts and renderer deleted)
- `tun-c36f` depends on `tun-5a52` (can start after core agent files deleted)
- `tun-be14` depends on `tun-c36f` (can start after constants cleaned up)

---

## Execution Notes

This is a straightforward deletion task following the same pattern as the planning feature deletion. The key insight is that the feature is **already disabled**, so there's no risk of breaking production. Work through milestones sequentially to catch import errors early.
