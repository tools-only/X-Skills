# Research – Codebase Bottleneck Analysis
**Date:** 2026-01-31
**Phase:** Research

## Executive Summary

This document maps the critical performance bottlenecks and architectural chokepoints in the tunacode codebase. The analysis reveals **10 high-priority bottleneck locations** spanning synchronous I/O in async contexts, shared mutable state without synchronization, O(n²) algorithms, and threading/async incompatibility issues.

---

## Structure Overview

```
src/tunacode/
├── ui/              # 197 Python files total
│   ├── app.py       # TextualReplApp - main event loop
│   ├── widgets/     # 10 UI components
│   └── renderers/   # 15+ output renderers
├── core/
│   ├── agents/
│   │   ├── main.py                    # RequestOrchestrator
│   │   ├── agent_components/
│   │   │   ├── streaming.py           # Streaming handler
│   │   │   ├── tool_executor.py       # Parallel execution
│   │   │   ├── orchestrator.py        # Node processing
│   │   │   └── response_state.py      # State with RLock
│   │   └── resume/
│   │       ├── sanitize.py            # Message cleanup
│   │       └── prune.py               # Token pruning
│   ├── state.py                       # StateManager singleton
│   └── types/tool_registry.py         # Mutable registry
└── tools/           # Tool implementations
```

---

## Critical Bottlenecks (Prioritized)

### 1. Synchronous File I/O in Async Context
**Location:** `src/tunacode/core/state.py:234-237`
```python
# save_session() - blocking call in async flow
with open(session_file, "w") as f:
    json.dump(session_data, f, indent=2)
```
**Impact:** Event loop stall during session persistence
**Pattern:** Called from `ui/app.py:250` in `_process_request()` finally block

**Related:**
- `load_session()` at `state.py:257-302` - similar blocking I/O
- `agent_config.py:174-175` - blocking file reads at startup

---

### 2. Threading Lock in Async Context (RLock)
**Location:** `src/tunacode/core/agents/agent_components/response_state.py:26`
```python
_lock: threading.RLock = field(default_factory=threading.RLock)
```
**Impact:** Blocks event loop when acquired from async contexts
**Usage:** All property accessors acquire lock (lines 52-59)

**Related:**
- `state_transition.py:49` - threading.RLock in StateMachine
- `core/logging/manager.py:22-35` - LogManager threading locks

---

### 3. Unbounded String Accumulation (O(n²) Memory)
**Location:** `src/tunacode/core/agents/agent_components/streaming.py:349-350`
```python
# Per-node debug accumulators - unbounded growth
state_manager.session._debug_raw_stream_accum += delta_text
```
**Impact:** Quadratic memory growth during large streaming responses
**Pattern:** String concatenation in loop (line 349-350, 153-154)

---

### 4. Shared Mutable State Without Synchronization
**Location:** `src/tunacode/core/types/tool_registry.py:33-40`
```python
@dataclass(slots=True)
class ToolCallRegistry:
    _calls: dict[ToolCallId, CanonicalToolCall] = field(default_factory=dict)
    _order: list[ToolCallId] = field(default_factory=list)
```
**Impact:** Race conditions during parallel tool execution
**Accessed by:** Main loop, tool executor, cleanup operations concurrently

---

### 5. In-Place List Mutation During Iteration
**Location:** `src/tunacode/core/agents/resume/sanitize.py:315`
```python
messages[:] = remaining_messages  # In-place mutation
```
**Impact:** Race conditions if agent holds references during cleanup
**Called from:** `main.py:461` during abort handling

**Related:**
- `run_cleanup_loop()` at `sanitize.py:458-502` - up to 10 iterations over messages
- O(n²) worst case from repeated list mutations

---

### 6. Non-Atomic Counter Increment (Race Condition)
**Location:** `src/tunacode/core/agents/agent_components/orchestrator/tool_dispatcher.py:298-299`
```python
batch_id = runtime.batch_counter + 1  # Read
runtime.batch_counter = batch_id       # Write (non-atomic)
```
**Impact:** ID collisions between concurrent tool batches
**TOCTOU Race:** State could change between read and write

---

### 7. Deep Attribute Access Chains (Violation of 2-Dot Rule)
**Location:** `src/tunacode/core/agents/main.py:96-100`
```python
session = self.state_manager.session
conversation = session.conversation
runtime = session.runtime
tool_registry = runtime.tool_registry  # 4 levels of nesting
```
**Impact:** High coupling, difficult refactoring, brittle code
**Pattern:** Violates "Law of Demeter" - reaches through 4+ object layers

**Related:**
- `ui/app.py:376-387` - UI layer directly accesses `session.conversation.total_tokens`

---

### 8. Sequential Batch Processing (Under-Utilized Parallelism)
**Location:** `src/tunacode/core/agents/agent_components/tool_executor.py:119-129`
```python
# Batches execute sequentially - batch 2 waits for batch 1
for i in range(0, len(tool_calls), max_parallel):
    batch = tool_calls[i:i + max_parallel]
    batch_tasks = [execute_with_retry(part, node) for part, node in batch]
    batch_results = await asyncio.gather(*batch_tasks)
```
**Impact:** If batch 1 has 4 slow tools, batch 2 waits despite available slots

---

### 9. Synchronous Token Counting (CPU Bottleneck)
**Location:** `src/tunacode/utils/messaging/token_counter.py`
**Called from:** `state.py:62-75` - `update_token_count()`
```python
for msg in messages:
    content = get_content(msg)
    token_count = estimate_tokens(content, model_name)  # Synchronous, O(n)
```
**Impact:** O(n) scan on every message, no caching per message

---

### 10. Debug Logging Overhead in Hot Path
**Location:** `src/tunacode/core/agents/agent_components/streaming.py:211-267`
```python
async for event in request_stream:
    debug_event_count += 1
    if debug_event_count <= DEBUG_STREAM_EVENT_LOG_LIMIT:
        # Complex reflection and logging executes even when debug_mode=False
```
**Impact:** Debug logic executes unconditionally (check is inside loop)

---

## Dependency Flow Violations

### UI → Core Boundary Violations
**Location:** `src/tunacode/ui/app.py:376-387`
```python
def _update_resource_bar(self) -> None:
    session = self.state_manager.session
    conversation = session.conversation
    context_tokens = conversation.total_tokens  # Deep reach into core state
```
**Issue:** UI layer reaches 3 levels deep into core state structures
**Correct pattern:** UI should use protocols, not direct state access

---

## Critical Request Flow (With Bottlenecks Marked)

```
User Input → Editor.action_submit() → widgets/editor.py:115-127
    └── EditorSubmitRequested message → ui/app.py:268-274
        └── request_queue.put() → ui/app.py:274
            └── _request_worker() → ui/app.py:177-186
                └── _process_request() → ui/app.py:188-251
                    └── process_request() → core/agents/main.py:553
                        └── RequestOrchestrator.run() → main.py:221-239
                            └── _run_impl() → main.py:241-482
                                └── agent.iter() → main.py:355
                                    └── stream_model_request_node() → streaming.py:122-384
                                        └── [BOTTLENECK: Unbounded string accumulation]
                                    └── process_node() → orchestrator.py:178-292
                                        └── dispatch_tools() → tool_dispatcher.py:219-352
                                            └── [BOTTLENECK: Non-atomic counter]
                                            └── execute_tools_parallel() → tool_executor.py:53-140
                                                └── [BOTTLENECK: Sequential batch processing]
                                                └── [BOTTLENECK: Shared registry access]
                        └── save_session() → state.py:234-237
                            └── [BOTTLENECK: Synchronous file I/O]
```

---

## Symbol Index – Key Functions/Classes

| Symbol | Location | Purpose | Risk |
|--------|----------|---------|------|
| `RequestOrchestrator` | `main.py:141` | Main request coordination | Single point of failure |
| `process_node()` | `orchestrator.py:178` | Per-node processing | 8 parameters, high coupling |
| `execute_tools_parallel()` | `tool_executor.py:53` | Parallel tool execution | Sequential batches |
| `ToolCallRegistry` | `tool_registry.py:33` | Tool call tracking | No synchronization |
| `ResponseState` | `response_state.py:26` | State management | RLock in async |
| `save_session()` | `state.py:234` | Session persistence | Blocking I/O |
| `run_cleanup_loop()` | `sanitize.py:458` | Message cleanup | O(n²) iterations |
| `prune_old_tool_outputs()` | `prune.py:160` | Token pruning | O(n) on large histories |
| `update_token_count()` | `state.py:62` | Token counting | No caching |
| `TextualReplApp` | `ui/app.py:79` | Main UI | Queue serializes requests |

---

## Patterns Found

### Anti-Pattern 1: Threading Locks in Async Code
**Locations:**
- `response_state.py:26` - threading.RLock
- `state_transition.py:49` - threading.RLock
- `logging/manager.py:22-35` - threading.RLock

**Fix:** Use `asyncio.Lock` or `asyncio.Semaphore` for async contexts

### Anti-Pattern 2: In-Place Mutation of Shared State
**Locations:**
- `sanitize.py:315` - `messages[:] = remaining_messages`
- `main.py:267-269` - `session_messages` mutated in-place

**Fix:** Create new lists, swap references atomically

### Anti-Pattern 3: Synchronous I/O in Async Functions
**Locations:**
- `state.py:234-237` - File write in `save_session()`
- `state.py:257-302` - File read in `load_session()`
- `agent_config.py:174-175` - File reads at startup

**Fix:** Wrap in `asyncio.to_thread()` or use `aiofiles`

### Anti-Pattern 4: Unbounded Accumulators
**Locations:**
- `streaming.py:349-350` - `_debug_raw_stream_accum`
- `streaming.py:153-154` - `_debug_events`

**Fix:** Limit size, use StringIO for accumulation, clear periodically

### Anti-Pattern 5: Deep Attribute Chains (4+ dots)
**Locations:**
- `main.py:96-100` - 4 levels to reach `tool_registry`
- `ui/app.py:376-387` - 3 levels to reach `total_tokens`

**Fix:** Extract needed values, pass as parameters, use protocols

---

## Contention Points

1. **StateManager.session** - Accessed by UI, core, and tools without coordination
2. **ToolCallRegistry** - Mutated during parallel tool execution without locks
3. **conversation.messages** - Mutated during cleanup while agent may hold references
4. **runtime.batch_counter** - Non-atomic increment under concurrent access
5. **LogManager singleton** - Threading locks don't yield to asyncio event loop

---

## Appendix: File Statistics

- **Total Python files:** 197
- **Core agent files:** ~25
- **UI files:** ~35
- **Tool files:** ~40
- **Test files:** Located in `tests/` directory

---

## Research Sources

- codebase-locator agent (0ab75708)
- codebase-analyzer agent (ade05d3b)
- context-synthesis agent (03730480)
- structure-map.sh output
- ast-scan.sh functions output
