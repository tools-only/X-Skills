# Research - Issue #311: core → types Layer Violation

**Date:** 2026-01-27
**Owner:** claude-agent
**Phase:** Research
**Issue:** https://github.com/alchemiststudiosDOTai/tunacode/issues/311

## Goal

Understand the scope and implications of the `core → types` layer violation reported in issue #311, and determine the correct remediation approach.

## Findings

### The Reported Problem

Issue #311 proposes a strict 5-layer hierarchy:

```
ui (Layer 0) → core (Layer 1) → tools/indexing/lsp (Layer 2) → utils (Layer 3) → types/configuration (Layer 4)
```

Under this model, `core` should NOT directly import from `types` (skips Layer 2 and 3).

### Current Architecture (DEPENDENCY_MAP.md)

The frozen baseline (2026-01-26) shows a different model:

```
ui          → outer layer (TUI)
core        → business logic
tools       → agent tools
indexing    → code indexing infrastructure
lsp         → language server protocol
─────────────────────────────────
utils       ┐
types       │ utils-level (importable from anywhere)
configuration│
constants   ┘
```

**Key difference:** Types is classified as "utils-level" - a horizontal layer that ANY vertical layer can import from. This is consistent with CLAUDE.md Gate 2.

### Architectural Contradiction

| Document | Position | Status |
|----------|----------|--------|
| Issue #311 | types is Layer 4, core cannot import | Proposed change |
| DEPENDENCY_MAP.md | types is utils-level, importable anywhere | Frozen baseline |
| CLAUDE.md Gate 2 | types/ is utils-level module | Authoritative |

The frozen baseline marks `core → types: 26 imports` as "✅ valid".

### Import Inventory (17 files, 26 imports)

#### State Management (5 imports)
| File | Types Imported |
|------|----------------|
| `core/state.py:17` | ConversationState, InputSessions, ModelName, RuntimeState, SessionId, UserConfig |
| `core/state.py:27` | UsageMetrics |
| `core/user_configuration.py:5-6` | UserConfig, StateManagerProtocol |
| `core/lsp_status.py:3` | UserConfig |

#### Core Types Re-export (2 imports)
| File | Types Imported |
|------|----------------|
| `core/types.py:5` | ModelName, ToolArgs, ToolCallback, ToolName, ToolResultCallback, ToolStartCallback |
| `core/types.py:14` | UsageMetrics |

#### Configuration (1 import)
| File | Types Imported |
|------|----------------|
| `core/configuration.py:31` | ModelPricing |

#### Logging (1 import)
| File | Types Imported |
|------|----------------|
| `core/logging/manager.py:10` | StateManagerProtocol |

#### Agent Core (1 import)
| File | Types Imported |
|------|----------------|
| `core/agents/main.py:22` | AgentRun, ModelName, NoticeCallback, StateManagerProtocol, StreamingCallback, ToolCallback, ToolResultCallback, ToolStartCallback |

#### Agent Components (8 imports)
| File | Types Imported |
|------|----------------|
| `core/agents/agent_components/agent_config.py:23` | ModelName, PydanticAgent, SessionStateProtocol |
| `core/agents/agent_components/agent_helpers.py:5` | CanonicalToolCall |
| `core/agents/agent_components/response_state.py:6` | AgentState |
| `core/agents/agent_components/state_transition.py:7` | AgentState |
| `core/agents/agent_components/streaming.py:16` | StreamingCallback |

#### Orchestrator (4 imports)
| File | Types Imported |
|------|----------------|
| `core/agents/agent_components/orchestrator/orchestrator.py:5-6` | AgentState, StreamingCallback, ToolCallback, ToolResultCallback, ToolStartCallback |
| `core/agents/agent_components/orchestrator/tool_dispatcher.py:12-13` | AgentState, ToolArgs, ToolCallId, ToolCallback, ToolResultCallback, ToolStartCallback |
| `core/agents/agent_components/orchestrator/usage_tracker.py:6-7` | UsageMetrics, normalize_request_usage |

#### Resume/Sanitization (4 imports)
| File | Types Imported |
|------|----------------|
| `core/agents/resume/sanitize.py:18-20` | ToolCallId, CanonicalMessage, MessageRole, SystemPromptPart, ToolCallRegistry |
| `core/agents/resume/sanitize_debug.py:10` | ToolCallId |

### Types Layer Analysis

The `src/tunacode/types/` directory contains 9 files with types in three categories:

#### Category A: Valid Utils-Level (Should Stay)

These types are used across multiple layers and are correctly placed:

| File | Types | Used By |
|------|-------|---------|
| `base.py` | 42 primitive aliases (ModelName, ToolArgs, etc.) | core, tools, utils |
| `pydantic_ai.py` | PydanticAgent, MessageHistory, NormalizedUsage | core, utils |
| `canonical.py` | CanonicalMessage, UsageMetrics, MessageRole, etc. | core, utils/messaging/adapter |
| `callbacks.py` | ToolCallback, StreamingCallback, etc. | core (dependency injection) |

#### Category B: Core-Specific (Candidates for Move)

These types are ONLY used by core and should not be in utils-level:

| File | Types | Reason |
|------|-------|--------|
| `dataclasses.py` | ResponseState | Only used in core/agents/agent_components/response_state.py |
| `dataclasses.py` | AgentState | Only used in core agent state transitions |
| `state.py` | SessionStateProtocol | Protocol for core's SessionState |
| `state.py` | StateManagerProtocol | Protocol for core's StateManager |
| `state_structures.py` | ConversationState, TaskState, RuntimeState, UsageState | Only compose core's SessionState |
| `tool_registry.py` | ToolCallRegistry | Stateful business logic, only used in RuntimeState and sanitize.py |

#### Category C: Dead Code (Delete)

| File | Types | Reason |
|------|-------|--------|
| `dataclasses.py` | CommandContext | Never imported anywhere |
| `callbacks.py` | ProcessRequestCallback | Only referenced in dead CommandContext |

### Key Patterns Found

**Pattern 1: Core Facade** (commit 10b48c30)
- When UI needed LSP functionality, a facade was created in core
- Pattern: `ui → core/facade.py → tools/implementation.py → lsp/`

**Pattern 2: Move to Utils** (commit be735754)
- Shared constants moved from tools to `utils/system/ignore_patterns.py`
- Both tools and indexing import from there

**Pattern 3: Central Re-export** (current)
- `types/__init__.py` re-exports everything, hiding internal structure
- `core/types.py` already re-exports a subset for core's use

## Knowledge Gaps

1. **Architectural Intent Unclear**: Is issue #311's strict layering the desired target state, or is the current "types as utils-level" the correct model?

2. **Migration Strategy**: If we adopt strict layering, do we:
   - Move core-specific types to `core/types/`?
   - Create re-exports in utils layer?
   - Use dependency injection?

3. **Backward Compatibility**: How do we handle the 26 existing imports without breaking changes?

## Recommendations

### Option A: Keep Current Architecture (Minimal Change)

If types/ remains "utils-level":
- Close issue #311 as "works as designed"
- Delete dead code (CommandContext, ProcessRequestCallback)
- Update documentation to clarify the model

### Option B: Strict Layer Enforcement (Issue #311's Proposal)

If we adopt strict layering:

**Phase 1: Move Core-Specific Types**
1. Create `src/tunacode/core/types/` directory
2. Move: ResponseState, AgentState, StateManagerProtocol, SessionStateProtocol
3. Move: ConversationState, TaskState, RuntimeState, UsageState
4. Move: ToolCallRegistry
5. Update 17 files in core to import from `core/types/`

**Phase 2: Cleanup**
1. Delete dead code: CommandContext, ProcessRequestCallback
2. Update `types/__init__.py` to only export utils-level types

**Phase 3: Verify**
1. Run grimp to confirm zero violations
2. Update DEPENDENCY_MAP.md

### My Recommendation

**Option B is the correct path forward**, but with nuance:

The current architecture has types/ as utils-level because it evolved organically. Issue #311 correctly identifies that core-specific types (protocols, state structures) should live in core, not in a shared layer.

However, the issue's proposed fix (re-exports through utils/tools) adds indirection. The cleaner fix is:

1. **Move core-specific types to `core/types/`** - types that describe core's internal structure belong in core
2. **Keep shared types in `types/`** - primitives, callbacks, canonical message types remain utils-level
3. **Result**: core imports its own types directly (no indirection), and shared types from utils-level (valid)

This maintains the "utils-level is importable anywhere" principle while moving core-specific implementation details out of the shared layer.

## References

- Issue: https://github.com/alchemiststudiosDOTai/tunacode/issues/311
- `docs/architecture/DEPENDENCY_MAP.md` - Current baseline
- `CLAUDE.md` Gate 2 - Dependency direction rules
- `src/tunacode/types/` - Types layer (9 files)
- `src/tunacode/core/types.py` - Existing core re-export pattern
