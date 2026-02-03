# Research – Recent Major Changes (Last 3 Days)
**Date:** 2026-01-28
**Owner:** Agent
**Phase:** Research

## Goal
Document all major changes from the last 3 days of development (Jan 25-28, 2026), including merged PRs, architectural refactors, and new features.

## Summary

The last 3 days saw **4 major PRs** merged with significant architectural improvements:

| PR | Title | Lines Changed | Key Impact |
|----|-------|---------------|------------|
| #319 | ChatContainer with insertion tracking | +1597 / -392 | Fixes tool panel race condition |
| #318 | Delete indexing system, fix LSP coupling | +228 / -1059 | Removes ~850 lines of dead code |
| #317 | Resolve layer violations (#313) | +1294 / -152 | Eliminates 24 core→utils violations |
| #316 | Core types layering | +587 / -various | Separates core/types from shared types |

---

## PR #319: ChatContainer with Insertion Tracking

**Merged:** 2026-01-28 00:06 UTC
**Author:** tunahorse

### Problem
Tool panels appeared at wrong position after stream cancellation or completion. Race condition: tool results arrive asynchronously via Textual's message queue, but `_current_stream` reference was cleared immediately on cancel/end.

### Solution
New `ChatContainer` widget (`src/tunacode/ui/widgets/chat.py`) with **insertion anchor tracking**:

1. **`start_stream()`** - Begins tracking, clears stale anchor
2. **`insert_before_stream()`** - Tool panels insert before streaming widget OR anchor
3. **`end_stream()`** - Captures finalized message as anchor for late panels
4. **`cancel_stream()`** - Preserves insertion context on cancel

### Key Files
- `src/tunacode/ui/widgets/chat.py` (217 lines) - New ChatContainer implementation
- `src/tunacode/ui/app.py` - Replaces RichLog with ChatContainer
- `src/tunacode/ui/welcome.py` - WriteableLog protocol for compatibility

### State Machine Pattern
| State | Condition | Insertion Behavior |
|-------|-----------|-------------------|
| **Streaming** | `_current_stream is not None` | Insert before streaming widget |
| **Finalized** | `_insertion_anchor is not None` | Insert before finalized message |
| **Idle** | Both `None` | Append to end |

---

## PR #318: Delete Indexing System and Fix LSP Lateral Coupling

**Merged:** 2026-01-27 19:49 UTC
**Author:** larock22

### Problem
1. Indexing system (~620 lines) was premature optimization with negligible benefit
2. Three lateral coupling violations between Layer 2 peers (tools, indexing, lsp)

### Solution

#### Deleted (~850 lines total)
- `src/tunacode/indexing/` - Entire module (code_index.py, constants.py)
- `src/tunacode/core/indexing_service.py` - Orchestration facade
- `src/tunacode/tools/lsp_status.py` - Merged into core
- LSP logic in `tools/decorators.py` - Removed orchestration from decorator

#### Added (~106 lines)
- `src/tunacode/tools/lsp/diagnostics.py` (66 lines) - Centralized `maybe_prepend_lsp_diagnostics()` helper
- Enhanced `src/tunacode/core/lsp_status.py` with `LspServerInfo` dataclass

### Lateral Coupling Fixes
| Violation | Fix |
|-----------|-----|
| tools/glob.py → indexing | Deleted indexing entirely |
| tools/decorators.py → lsp | Tools now call diagnostics helper directly |
| tools/lsp_status.py → lsp | Merged into core/lsp_status.py |

### New LSP Diagnostics Flow
```
write_file/update_file
  └─> maybe_prepend_lsp_diagnostics(result, filepath)
      └─> get_diagnostics() with timeout
      └─> format_diagnostics() as XML
      └─> Prepend to tool result
```

---

## PR #317: Resolve Layer Violations (fixes #313)

**Merged:** 2026-01-27 19:08 UTC
**Author:** larock22

### Problem
24 direct imports from `core` to `utils` bypassing the tools layer.

### Solution: Module Reorganization

#### Moved to `configuration/` (Foundation Layer)
- `utils/config/user_configuration.py` → `configuration/user_config.py`
- `utils/limits.py` → `configuration/limits.py`
- `utils/system/paths.py` → `configuration/paths.py`
- `utils/system/ignore_patterns.py` → `configuration/ignore_patterns.py`

#### Moved to `tools/` (Layer 2)
- `utils/parsing/` → `tools/parsing/` (command_parser, tool_parser, json_utils, retry)

#### Moved to `infrastructure/`
- `utils/ui/file_filter.py` → `infrastructure/file_filter.py` (standalone, zero internal deps)

#### Created Facade
- `tools/messaging/__init__.py` - Re-exports from `utils/messaging/`

### Dependency Impact
| From → To | Before | After |
|-----------|--------|-------|
| core → utils | 24 imports | **0 imports** |
| tools → utils | 4 imports | 2 imports |

---

## PR #316: Core Types Layering

**Merged:** 2026-01-27 (earlier)
**Author:** larock22

### Problem
Core-specific state types lived in shared `types/` foundation layer, creating namespace pollution.

### Solution

#### Created `core/types/` Package
Exports core-only types:
- `AgentState` enum (USER_INPUT, ASSISTANT, TOOL_EXECUTION, RESPONSE)
- `ResponseState` dataclass
- `SessionStateProtocol`, `StateManagerProtocol`
- `ConversationState`, `RuntimeState`, `TaskState`, `UsageState`

#### Created `core/shared_types.py` Facade
UI-accessible types only:
- `ModelName`, `ToolArgs`, `ToolCallback`, `ToolName`
- `UserConfig`, `UsageMetrics`

### Benefit
UI sees only what core chooses to expose via facade. Core internals protected.

---

## Other Notable Commits

### LSP Server Status in ResourceBar
- **Commit:** `a75d270b`
- **File:** `src/tunacode/ui/widgets/resource_bar.py`
- Shows LSP server status (running/stopped) in UI status bar

### File Filter Fuzzy Matching Fix
- **Commit:** `087025ce`
- **File:** `src/tunacode/infrastructure/file_filter.py`
- Fixed fuzzy matching in file autocomplete

### Dependency Layer Visualization
- **Commits:** `74aa87d4`, `601a7c1c`
- **Files:** `scripts/grimp_layers_report.py`, `docs/architecture/DEPENDENCY_LAYERS.png`
- Automated PNG generation from layer report

### Mypy Batch Fixes
- **Commit:** `a77fd94c`
- Resolved multiple mypy type errors across configuration, core, and tools

---

## Architectural State After Changes

### Layer Hierarchy (Enforced)
```
ui (Layer 0)
    ↓
core (Layer 1)
    ↓
infrastructure (plugin)
    ↓
tools (Layer 2)
    ↓
utils (Layer 3)
    ↓
[Foundation: configuration, types, constants, exceptions]
```

### Key Patterns Established

1. **Foundation Layer Promotion** - Static data/pure functions → `configuration/`
2. **Layer 2 Re-export Facade** - `tools/messaging/` delegates to `utils/messaging/`
3. **Domain Type Extraction** - `core/types/` for core-only, `core/shared_types.py` for UI
4. **Infrastructure Plugin** - Standalone modules wrapped by core facades

### Verification
```bash
# Zero layer violations
grep -r "from tunacode.utils" src/tunacode/core/ | wc -l  # Returns 0
```

---

## Test Status
- All 512 existing tests pass
- Type checks clean (no new mypy errors)
- Linters pass

---

## References

### Research Documents
- `memory-bank/research/2026-01-27_11-53-37_issue-313-core-utils-layer-violation.md`
- `memory-bank/research/2026-01-27_insert-before-stream-investigation.md`
- `memory-bank/research/2026-01-27_17-45-00_chatcontainer_refactor.md`

### Architecture Documents
- `docs/architecture/DEPENDENCY_MAP.md`
- `docs/architecture/DEPENDENCY_LAYERS.md`
- `docs/lsp-diagnostics.md`

### Execution Logs
- `memory-bank/execute/2026-01-27_17-52-14_insert-before-stream-fix.md`
- `memory-bank/execute/2026-01-27_21-50-00_tickets.md`
