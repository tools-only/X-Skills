# Research – Tool Dispatcher Complexity (#358)
**Date:** 2026-02-05
**Owner:** Claude Code
**Phase:** Completed (see PR #366)
**Git Commit:** 9ed9739dda91c3f4f463c4c8f58687f75879afc7
**Completed:** 2026-02-05 via PR #366

> **Note:** This research was completed in PR #366. The tool_dispatcher was decomposed into 6 focused submodules. See the plan document and PR for implementation details.

## Goal
Map the current complexity issues in `tool_dispatcher.py` and identify refactoring patterns in the codebase for reducing complexity while maintaining behavior.

## Findings

### Current File Metrics
| Metric | Value |
|--------|-------|
| Total lines | 464 |
| Functions | 15 |
| Deferred imports | 5 |
| Magic numbers | 3 instances of `100` (`[:100]` debug preview truncation) + 1 instance of `1000` (ms conversion) |
| Debug mode conditionals | 6 `if` statements guarded by `debug_mode` |

### Complexity Hotspots

#### `dispatch_tools` function (lines 413-464)
- **Lines:** 52 total; 30 excluding docstring
- **Branches / decision points:** 4 `if` statements; ~7 decision points if counting `and`/boolean ops
- **Nesting levels:** 2 maximum (file max is 3 in `_collect_structured_tool_calls`)
- **Phases:** 4 sequential phases (collect structured → fallback → execute → log summary) + optional state-transition hook

#### `_collect_structured_tool_calls` (lines 220-267)
- **Complexity:** ~8 (loop + 4 conditionals + 2 debug branches)
- Dead parameters: `tool_callback` (line 224), `node` (line 222) - passed but never used

#### `_collect_fallback_tool_calls` (lines 275-341)
- **Complexity:** ~10 (loop + 5 conditionals + debug branch + assertions)
- Contains assertions: `assert isinstance(result, tuple)` and `assert isinstance(result, list)`

### Dead Code Identified

| Pattern | Location | Details |
|---------|----------|---------|
| Unused parameters | Lines 222, 224 | `_collect_structured_tool_calls` accepts `node` and `tool_callback` but never references them |
| Redundant local (derived output) | Lines 437, 446, 464 | `used_fallback` does not affect control flow in `dispatch_tools`, but it is part of `ToolDispatchResult` and asserted in tests (so not removable without an API change) |

### Magic Numbers

| Value | Lines | Context | Recommendation |
|-------|-------|---------|----------------|
| `100` | 257, 259, 306 | Debug preview string slicing (`[:100]`) | Extract to `DEBUG_PREVIEW_MAX_LENGTH` |
| `1000` | 461 | Seconds → milliseconds conversion | Extract to `MS_PER_SECOND` |

### Deferred Imports (5 total)

| Function | Imports | Location |
|----------|---------|----------|
| `normalize_tool_args` | `parse_args` | `tunacode.tools.parsing.command_parser` |
| `_ensure_normalized_tool_call_part` | `ToolCallPart` | `pydantic_ai.messages` |
| `_collect_fallback_tool_calls` | `ToolCallPart`, `has_potential_tool_call`, `parse_tool_calls_from_text` | `pydantic_ai.messages`, `tunacode.tools.parsing.tool_parser` |
| `_execute_tool_batch` | `execute_tools_parallel` | `..tool_executor` |

### Seven Responsibilities Mapped

| # | Responsibility | Line Range | Function(s) |
|---|----------------|------------|-------------|
| 1 | **Collect structured tool calls** | 220-267 | `_collect_structured_tool_calls` |
| 2 | **Collect fallback tool calls** | 275-341 | `_collect_fallback_tool_calls` |
| 3 | **Register/normalize tool names** | 60-65, 73-83, 192-212 | `_normalize_tool_name`, `_register_tool_call`, `_ensure_normalized_tool_call_part` |
| 4 | **Parse and normalize tool arguments** | 127-131, 134-155 | `normalize_tool_args`, `record_tool_call_args` |
| 5 | **Mark tool calls as running** | 85-94 | `_mark_tool_calls_running` |
| 6 | **Execute tools in parallel** | 349-380 | `_execute_tool_batch` |
| 7 | **Record tool failures** | 97-119 | `_record_tool_failure` |

Additional responsibilities:
- **Debug logging/diagnostics** - Scattered throughout (lines 246-265, 288, 300-327)
- **Tool name validation** - Lines 51-57, 60-65
- **Dispatch summary logging** - Lines 388-405
- **Tool call args retrieval** - Lines 158-166

### Related Files

**Files using tool_dispatcher.py:**
- `/home/tuna/tunacode/src/tunacode/core/agents/agent_components/orchestrator/orchestrator.py` - Imports `dispatch_tools`, `has_tool_calls`
- `/home/tuna/tunacode/src/tunacode/core/agents/agent_components/orchestrator/tool_returns.py` - Imports `consume_tool_call_args`

**Test files:**
- `/home/tuna/tunacode/tests/integration/tools/test_tool_dispatcher_coverage.py` - 10 test functions
- `/home/tuna/tunacode/tests/integration/core/test_tool_call_lifecycle.py` - 26 test functions
- Total: 36 test functions across both files

**Sibling orchestrator modules:**
- `orchestrator.py` - Main orchestrator
- `tool_returns.py` - Tool return processing
- `message_recorder.py` - Message recording
- `usage_tracker.py` - Token usage tracking
- `debug_format.py` - Debug formatting

## Key Patterns / Solutions Found

### 1. Node Orchestrator Pattern (orchestrator.py)
The main orchestrator delegates each step to focused modules:
1. Transition to ASSISTANT state
2. Process inbound tool returns (`emit_tool_returns`)
3. Record agent thought (`record_thought`)
4. Track token usage (`update_usage`)
5. Extract text content
6. **Dispatch outbound tool calls (`dispatch_tools`)**
7. Check for node result
8. Transition to RESPONSE state

This is the pattern tool_dispatcher.py should follow - further decomposition into focused modules.

### 2. 4-Phase Dispatch Architecture
```python
# Phase 1: Collect structured tool calls
# Phase 2: Fallback to text parsing if no structured calls
# Phase 3: Execute if we have a callback
# Phase 4: Log summary
```
These phases could be extracted into separate modules or classes.

### 3. Registry Pattern (tool_registry.py)
The `ToolCallRegistry` class at `/home/tuna/tunacode/src/tunacode/core/types/tool_registry.py` demonstrates:
- Immutable updates via `dataclasses.replace()`
- Order preservation for listing
- Lifecycle tracking: PENDING → RUNNING → COMPLETED/FAILED/CANCELLED

### 4. Tool Executor Pattern (tool_executor.py)
- Batched parallel execution respecting `TUNACODE_MAX_PARALLEL`
- Exponential backoff with jitter for retries
- Non-retryable error whitelist
- Failure callback notification

### 5. State Machine Pattern (state_transition.py)
- Thread-safe state transitions with explicit rules
- State machine ownership stays in orchestrator
- Transitions injected via callbacks (e.g., `response_state_transition`)

## Knowledge Gaps

1. **Test Coverage:** Are there edge cases in `_collect_fallback_tool_calls` that are not covered by tests?
2. **Debug Mode Usage:** The 6 `debug_mode`-guarded conditionals add branching complexity - are all necessary?
3. **Fallback Parsing:** Is the fallback text parsing (`_collect_fallback_tool_calls`) actively used, or legacy code?
4. **Import Coupling:** Why are imports deferred inside functions? (likely to avoid circular imports)
5. **State Transition Callback:** The `response_state_transition` callback is passed but not always called - clarify contract

## Refactoring Opportunities

### Immediate (Low Risk)
1. Extract magic number `100` → `DEBUG_PREVIEW_MAX_LENGTH`
2. Remove dead parameters `node` and `tool_callback` from `_collect_structured_tool_calls`
3. Add constant for `MS_PER_SECOND = 1000`

### Medium (Behavior Preserving)
1. **Extract Tool Collection Module:** Move `_collect_structured_tool_calls` and `_collect_fallback_tool_calls` to a new `tool_collection.py` module
2. **Extract Tool Registry Module:** Move `_register_tool_call`, `_mark_tool_calls_running`, `_record_tool_failure` to a registry operations module
3. **Extract Logging Module:** Move `_log_dispatch_summary` and debug logging to a separate module
4. **Extract Constants:** Move all constants at top of file to a separate file

### Architectural
1. **Create ToolDispatcher Class:** Instead of module-level functions, create a class to hold state and reduce parameter passing
2. **Simplify `dispatch_tools`:** The main function should orchestrate, not implement - delegate all work to sub-modules

## References

### GitHub Permalinks
- tool_dispatcher.py: https://github.com/tunacode/tunacode/blob/9ed9739dda91c3f4f463c4c8f58687f75879afc7/src/tunacode/core/agents/agent_components/orchestrator/tool_dispatcher.py
- orchestrator.py: https://github.com/tunacode/tunacode/blob/9ed9739dda91c3f4f463c4c8f58687f75879afc7/src/tunacode/core/agents/agent_components/orchestrator/orchestrator.py
- tool_registry.py: https://github.com/tunacode/tunacode/blob/9ed9739dda91c3f4f463c4c8f58687f75879afc7/src/tunacode/core/types/tool_registry.py
- tool_executor.py: https://github.com/tunacode/tunacode/blob/9ed9739dda91c3f4f463c4c8f58687f75879afc7/src/tunacode/core/agents/agent_components/tool_executor.py

### Test Files
- test_tool_dispatcher_coverage.py: https://github.com/tunacode/tunacode/blob/9ed9739dda91c3f4f463c4c8f58687f75879afc7/tests/integration/tools/test_tool_dispatcher_coverage.py
- test_tool_call_lifecycle.py: https://github.com/tunacode/tunacode/blob/9ed9739dda91c3f4f463c4c8f58687f75879afc7/tests/integration/core/test_tool_call_lifecycle.py

### Issue
- GitHub Issue #358: Tool dispatcher complexity and dead code
