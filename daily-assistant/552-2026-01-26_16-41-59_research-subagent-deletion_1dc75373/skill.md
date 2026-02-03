# Research - Research Sub-Agent Feature Deletion

**Date:** 2026-01-26
**Owner:** agent
**Phase:** Research

## Goal

Map all components of the `research_codebase` sub-agent feature for clean removal. The feature is already disabled in production (commented out in agent_config.py).

## Findings

### Critical Status: ALREADY DISABLED

The research tool is already commented out in production:

```python
# src/tunacode/core/agents/agent_components/agent_config.py:37-39
# TODO: Re-enable research_codebase subagent after fixing parameter name mismatch
# See: memory-bank/research/2026-01-25_limited-read-file-parameter-mismatch.md
# from tunacode.core.agents.delegation_tools import create_research_codebase_tool
```

**Impact**: No active users, no runtime impact, safe to delete.

---

## Files to DELETE (7 files)

### Core Implementation (2 files)

| File | Lines | Purpose |
|------|-------|---------|
| `src/tunacode/core/agents/delegation_tools.py` | 154 | Tool factory with `create_research_codebase_tool()`, `ResearchResult` TypedDict |
| `src/tunacode/core/agents/research_agent.py` | 249 | Agent factory with `create_research_agent()`, `ProgressTracker`, limited read_file wrapper |

### UI Renderer (1 file)

| File | Lines | Purpose |
|------|-------|---------|
| `src/tunacode/ui/renderers/tools/research.py` | 296 | `ResearchRenderer` class, `ResearchData` dataclass, 4-zone panel layout |

### Prompt Sections (4 files)

| File | Purpose |
|------|---------|
| `src/tunacode/prompts/research/sections/agent_role.xml` | Research agent role definition |
| `src/tunacode/prompts/research/sections/tool_use.xml` | Tool usage guidelines |
| `src/tunacode/prompts/research/sections/constraints.xml` | Research constraints |
| `src/tunacode/prompts/research/sections/output_format.xml` | Output format spec |

Also delete the parent directory: `src/tunacode/prompts/research/`

---

## Files to MODIFY (17 files)

### Core Layer (5 files)

| File | Line(s) | Change |
|------|---------|--------|
| `src/tunacode/core/agents/agent_components/agent_config.py` | 37-39, 423-428 | Remove commented import and tool registration |
| `src/tunacode/core/agents/agent_components/orchestrator/tool_dispatcher.py` | 23, 30 | Remove `RESEARCH_TOOL_NAME` and `TOOL_START_RESEARCH_LABEL` constants |
| `src/tunacode/core/agents/agent_components/agent_helpers.py` | 54 | Remove research case from `get_readable_tool_description()` |
| `src/tunacode/core/prompting/templates.py` | 48-60 | Remove `RESEARCH_TEMPLATE` definition |
| `src/tunacode/core/prompting/__init__.py` | 14, 26 | Remove `RESEARCH_TEMPLATE` import and export |

### Constants & Config (2 files)

| File | Line(s) | Change |
|------|---------|--------|
| `src/tunacode/constants.py` | 78, 88 | Remove `RESEARCH_CODEBASE` from `ToolName` enum and `READ_ONLY_TOOLS` |
| `src/tunacode/configuration/settings.py` | 32 | Remove `ToolName.RESEARCH_CODEBASE` from tool list |

### UI Layer (2 files)

| File | Line(s) | Change |
|------|---------|--------|
| `src/tunacode/ui/renderers/tools/__init__.py` | 29, 69 | Remove `render_research_codebase` import and `__all__` export |
| `src/tunacode/ui/renderers/panels.py` | 504, 518 | Remove research renderer from registry |

### Tests (2 files)

| File | Line(s) | Change |
|------|---------|--------|
| `tests/integration/tools/test_tool_dispatcher_coverage.py` | 30, 62-65, 236 | Remove `RESEARCH_TOOL_NAME` constant and test |
| `tests/integration/core/test_tool_call_lifecycle.py` | 75 | Remove research from filter logic |

### Scripts (1 file)

| File | Line(s) | Change |
|------|---------|--------|
| `scripts/preview_tool_panels.py` | 211, 264, 390-396 | Remove research preview builder |

### Documentation (5 files)

| File | Section | Change |
|------|---------|--------|
| `docs/codebase-map/modules/core-agents.md` | 102-154 | Remove delegation/research section |
| `docs/codebase-map/modules/core-prompting.md` | 34, 61 | Remove `RESEARCH_TEMPLATE` documentation |
| `docs/codebase-map/modules/ui-overview.md` | 90 | Remove research renderer reference |
| `docs/codebase-map/structure/02-core-directory.md` | 30, 60-64 | Remove research_agent.py mentions |
| `docs/codebase-map/structure/tree-structure.txt` | 63, 160 | Remove research files from tree |

---

## Key Patterns / Solutions Found

### Architecture Layers Affected

```
┌─────────────────────────────────────────────────┐
│  PROMPTS         research/sections/ (4 XML)     │  DELETE
├─────────────────────────────────────────────────┤
│  TEMPLATES       RESEARCH_TEMPLATE              │  REMOVE
├─────────────────────────────────────────────────┤
│  AGENTS          delegation_tools.py            │  DELETE
│                  research_agent.py              │  DELETE
├─────────────────────────────────────────────────┤
│  CONSTANTS       ToolName.RESEARCH_CODEBASE     │  REMOVE
│                  READ_ONLY_TOOLS                │  REMOVE
├─────────────────────────────────────────────────┤
│  DISPATCHER      RESEARCH_TOOL_NAME             │  REMOVE
├─────────────────────────────────────────────────┤
│  UI              research.py (renderer)         │  DELETE
│                  panels.py (registry)           │  MODIFY
├─────────────────────────────────────────────────┤
│  TESTS           dispatcher_coverage.py         │  MODIFY
│                  tool_call_lifecycle.py         │  MODIFY
└─────────────────────────────────────────────────┘
```

### Comparison to Planning Feature Deletion

| Layer | Planning (completed) | Research (this task) |
|-------|---------------------|----------------------|
| Tool | `present_plan.py` | `delegation_tools.py` + `research_agent.py` |
| UI | `plan_approval.py` | `research.py` (renderer) |
| Prompts | `present_plan_prompt.xml` | `prompts/research/sections/` (4 files) |
| Constants | `EXIT_PLAN_MODE_SENTINEL` | `RESEARCH_TOOL_NAME`, `MAX_RESEARCH_FILES` |
| Types | `PlanApprovalCallback` | `ResearchResult`, `ResearchCodeExample` |
| Status | Was active | Already disabled |

---

## Knowledge Gaps

1. **Circular import handling**: `research_agent.py:194` has lazy import from `agent_config.py` - verify no other modules depend on this pattern
2. **Test coverage**: Need to verify no additional tests import from `delegation_tools.py` or `research_agent.py`
3. **Script dependencies**: `preview_tool_panels.py` may have additional research references beyond line 390

---

## Deletion Strategy (5 Phases)

### Phase 1: Delete Standalone Files
- Delete prompt sections directory
- Delete research renderer

### Phase 2: Remove Agent Files
- Delete `delegation_tools.py`
- Delete `research_agent.py`

### Phase 3: Remove Constants & Config
- Remove from `ToolName` enum
- Remove from `READ_ONLY_TOOLS`
- Remove from settings
- Remove from dispatcher constants

### Phase 4: Remove Prompting Template
- Remove `RESEARCH_TEMPLATE` definition
- Remove exports from `__init__.py`

### Phase 5: Update Tests & Docs
- Update test files
- Update documentation
- Update preview script

---

## References

- `src/tunacode/core/agents/delegation_tools.py` - Main tool factory
- `src/tunacode/core/agents/research_agent.py` - Agent creation
- `src/tunacode/ui/renderers/tools/research.py` - UI renderer
- `memory-bank/research/2026-01-25_limited-read-file-parameter-mismatch.md` - Bug that disabled feature
- `memory-bank/research/2026-01-26_planning-feature-deletion.md` - Similar deletion pattern
