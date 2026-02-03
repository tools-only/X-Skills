# Research – Agent Bottlenecks Analysis
**Date:** 2026-02-02 20:47:34
**Owner:** Claude Code (Research)
**Phase:** Research
**Git Commit:** 4fa9eb15d09e1263d0bd443f4ac7469b52d9904b
**Last Updated:** 2026-02-02 20:47:34
**Last Updated By:** Claude Code (Research)
**Tags:** [performance, bottlenecks, agent, tools, streaming]

## Goal
Identify major bottlenecks that can slow down the agent, excluding the known streaming issue.

---

## Additional Search
- `grep -ri "bottleneck\|performance\|overhead" .claude/`

---

## Findings

### Relevant Files & Why They Matter

#### Critical Bottlenecks

| File | Lines | Issue | Severity |
|------|-------|-------|----------|
| `src/tunacode/core/agents/resume/sanitize.py` | 458-502 | Redundant message canonicalization in cleanup loop | **HIGH** |
| `src/tunacode/core/agents/resume/prune.py` | 161-238 | Token estimation overhead in backward scan | **HIGH** |
| `src/tunacode/tools/parsing/command_parser.py` | 32-40 | JSON parsing with up to 10 retries | **HIGH** |
| `src/tunacode/core/agents/agent_components/tool_executor.py` | 119-130 | Sequential batch processing | **MEDIUM** |
| `src/tunacode/core/agents/agent_components/orchestrator/tool_dispatcher.py` | 235-368 | dispatch_tools() McCabe 29 complexity | **MEDIUM** |
| `src/tunacode/core/agents/main.py` | 303, 342-380 | Redundant logging even when debug_mode=False | **MEDIUM** |

---

### 1. Redundant Message Canonicalization (CRITICAL)

**Location:** `src/tunacode/core/agents/resume/sanitize.py:458-502`

**Issue:** The cleanup loop runs up to 10 iterations (`MAX_CLEANUP_ITERATIONS`). Each iteration calls multiple functions that independently convert all messages to canonical format:

```python
# Line 472: find_dangling_tool_call_ids() converts to canonical
dangling_tool_call_ids = find_dangling_tool_call_ids(messages)

# Line 473: remove_dangling_tool_calls() - no conversion here
removed_dangling = remove_dangling_tool_calls(messages, dangling_tool_call_ids, ...)

# Line 483: remove_empty_responses() converts to canonical AGAIN
removed_empty_responses = remove_empty_responses(messages)

# Line 488: remove_consecutive_requests() converts to canonical AGAIN
removed_consecutive_requests = remove_consecutive_requests(messages)
```

**Performance Impact:**
- For a conversation with 1000 messages: up to 30,000 message conversions (10 iterations × 1000 messages × 3 conversions)
- Each conversion involves nested loops through message parts (`adapter.py:146-202`)
- Runs on EVERY request before the agent starts

**Fix:** Convert messages to canonical format ONCE at the start of `run_cleanup_loop()`, then pass the canonical list to child functions.

---

### 2. Token Estimation Overhead in Pruning (CRITICAL)

**Location:** `src/tunacode/core/agents/resume/prune.py:161-238`

**Issue:** The `prune_old_tool_outputs()` function performs token estimation on EVERY tool return part during a backward scan:

```python
# Lines 192-202
for msg_idx in range(len(messages) - 1, -1, -1):
    message = messages[msg_idx]
    if not hasattr(message, "parts"):
        continue
    parts = message.parts
    for part_idx in range(len(parts) - 1, -1, -1):
        part = parts[part_idx]
        if is_tool_return_part(part):
            tokens = estimate_part_tokens(part, model_name)  # Called for each part
```

**Performance Impact:**
- O(M * P * L) where M = messages, P = parts per message, L = content length
- For large conversations with many tool calls, this processes hundreds of parts
- Token estimation uses string length division (`estimate_tokens()`)

**Fix:** Use cached token counts per part, or switch to byte-length approximation instead of character division.

---

### 3. JSON Parsing Retry Overhead (HIGH)

**Location:** `src/tunacode/tools/parsing/command_parser.py:32-40`

**Issue:** Tool argument parsing uses exponential backoff with up to 10 retries:

```python
# Lines 32-40
async def parse_args(args_str: str) -> dict[str, Any]:
    return await retry_json_parse_async(
        args_str,
        max_retries=10,      # 10 retries!
        base_delay=0.1,      # Starts at 100ms
        max_delay=5.0,       # Max 5 second delay
    )
```

**Performance Impact:**
- Failed JSON parse can cause delays up to 5+ seconds per tool call
- Multiple tools with JSON errors compound the delay
- Called sequentially before parallel tool execution begins

**Fix:** Consider reducing max_retries for tool args, or using a fast-path for common error patterns.

---

### 4. Sequential Batch Processing (MEDIUM)

**Location:** `src/tunacode/core/agents/agent_components/tool_executor.py:119-130`

**Issue:** When tool count exceeds `max_parallel`, batches execute sequentially:

```python
# Lines 120-130
if len(tool_calls) > max_parallel:
    results: list[None] = []
    for i in range(0, len(tool_calls), max_parallel):
        batch = tool_calls[i : i + max_parallel]
        batch_tasks = [execute_with_retry(part, node) for part, node in batch]
        batch_results = await asyncio.gather(*batch_tasks, return_exceptions=True)
        # Next batch waits for this batch to complete
```

**Performance Impact:**
- With 20 tool calls and max_parallel=5: 4 sequential batches
- If batch 1 completes quickly but batch 2 has slow tools, batch 3 waits unnecessarily
- Under-utilizes available parallelism

**Fix:** Use a task queue with `max_parallel` concurrent workers instead of fixed batches.

---

### 5. High McCabe Complexity in dispatch_tools (MEDIUM)

**Location:** `src/tunacode/core/agents/agent_components/orchestrator/tool_dispatcher.py:235-368`

**Issue:** The `dispatch_tools()` function has McCabe complexity of 29 (threshold is 25):

**Performance Impact:**
- Complex control flow makes optimization difficult
- Multiple nested conditionals and loops
- Hard to maintain and reason about performance characteristics

**Fix:** Refactor into smaller, focused functions with single responsibilities.

---

### 6. Redundant Logging Overhead (MEDIUM)

**Location:** `src/tunacode/core/agents/main.py:303, 342-380`

**Issue:** `_log_history_state()` performs detailed message inspection on EVERY request, even when `debug_mode` is False:

```python
# Lines 363-379: Called every request
message_sample = session_messages[-DEBUG_HISTORY_SAMPLE_SIZE:]
for idx, msg in enumerate(message_sample):
    msg_kind = getattr(msg, "kind", "?")
    msg_parts = getattr(msg, "parts", [])
    parts_summary = []
    parts_to_log = msg_parts[:DEBUG_HISTORY_PARTS_LIMIT]
    for part in parts_to_log:
        part_kind = getattr(part, "part_kind", "?")
        part_content = getattr(part, "content", None)
        # ... more processing
```

**Performance Impact:**
- Loops through last 3 messages
- For each message, loops through up to 5 parts
- Multiple `getattr()` calls per part
- Runs even when `debug_mode` is False

**Fix:** Only call `_log_history_state()` when `debug_mode` is True.

---

### 7. Unnecessary Threading Locks (LOW)

**Location:** `src/tunacode/core/agents/agent_components/response_state.py:49-111`

**Issue:** Every boolean property access acquires a reentrant lock:

```python
# Lines 50-59
@property
def has_user_response(self) -> bool:
    with self._lock:  # Lock in single-threaded asyncio
        return self._has_user_response
```

**Performance Impact:**
- Lock acquisition overhead on every property access
- asyncio is single-threaded - locks are unnecessary
- Locks only add overhead without providing benefit

**Fix:** Remove threading locks since the code runs in a single asyncio event loop.

---

### 8. Tool Registry List Operations (LOW)

**Location:** `src/tunacode/core/types/tool_registry.py:144-158`

**Issue:** `remove_many()` calls `remove()` in a loop, each doing an O(N) list search:

```python
# Lines 152-158
def remove_many(self, tool_call_ids: Iterable[ToolCallId]) -> int:
    removed_count = 0
    for tool_call_id in tool_call_ids:
        if self.remove(tool_call_id):  # O(N) per call
            removed_count += 1
    return removed_count
```

**Performance Impact:**
- O(N * M) where N = registry size, M = items to remove
- Called during cleanup with potentially many dangling tool call IDs

**Fix:** Build a set of IDs to remove, then rebuild `_order` list in one pass.

---

### 9. Unused/Dead Code Issues

#### 9.1 should_compact() Scans All Messages (Unused)

**Location:** `src/tunacore/core/agents/resume/summary.py:97-125`

**Issue:** The `should_compact()` function scans ALL messages to estimate total tokens, but appears to be unused in the main flow.

**Evidence:**
- Lines 113-123: Nested loops through all messages and parts
- Calls `estimate_tokens()` on every text part
- Exported but never called in `main.py` or `app.py`
- 503 lines in `summary.py` that may be dead code

**Fix:** Either wire `should_compact()` into the main flow OR mark as experimental and document as unused.

---

#### 9.2 Unused total_tokens Field

**Location:** `src/tunacode/core/types/state_structures.py:30`

**Issue:** `ConversationState.total_tokens` is maintained but never used after the recent refactor.

**Evidence:**
- Delta card `2026-02-02-remove-heuristic-token-counting.md` states: "Removed `update_token_count` calls"
- `total_tokens` field still exists in the dataclass
- UI now uses `session.usage.session_total_usage.total_tokens` instead

**Fix:** Deprecate `ConversationState.total_tokens` since it's no longer the authoritative source.

---

### 10. Fallback Tool Parser Overhead (LOW-MEDIUM)

**Location:** `src/tunacode/tools/parsing/tool_parser.py:264-323`

**Issue:** Four sequential regex-based parsing strategies are tried:

```python
# Lines 244-249
PARSING_STRATEGIES = [
    parse_qwen2_xml,
    parse_hermes_style,
    parse_code_fence,
    parse_raw_json,
]
```

**Performance Impact:**
- Each strategy runs regex matching on the full text
- 4 strategies × O(text length) per failed parse
- Only runs when primary tool call extraction fails

**Fix:** Consider caching parsed results, or use a pre-filter to skip unlikely strategies.

---

## Key Patterns / Solutions Found

### Performance Anti-Patterns Identified

1. **Repeated full-history scans** - Multiple functions independently iterate through all messages
2. **Redundant format conversions** - Canonical conversion happens multiple times per request
3. **Unnecessary synchronization** - Threading locks in single-threaded asyncio context
4. **Sequential batching** - Under-utilizes available parallelism
5. **High complexity functions** - McCabe score > 25 makes optimization difficult

### Optimization Opportunities

| Priority | Fix | Expected Impact |
|----------|-----|-----------------|
| 1 | Cache canonical conversion in cleanup loop | 10x reduction in message processing |
| 2 | Optimize prune token counting | 30-50% faster pruning |
| 3 | Move _log_history_state inside debug_mode check | Eliminate unnecessary inspection |
| 4 | Use task queue for tool execution | Better parallelism utilization |
| 5 | Remove threading locks | Eliminate lock overhead |

---

## Knowledge Gaps

1. **Actual performance measurements** - No profiling data to quantify actual time spent in each bottleneck
2. **Typical conversation size** - Don't know average message/part counts in real usage
3. **Tool call frequency** - Don't know typical number of tool calls per request
4. **JSON error rate** - Don't know how often JSON parsing fails and triggers retries

---

## Complexity Hotspots

From previous research (`memory-bank/research/2026-02-01_complexity-hotspots.md`):

| Function | File | McCabe | Threshold | Status |
|----------|------|--------|-----------|--------|
| `dispatch_tools` | `tool_dispatcher.py` | 29 | 25 | Exceeds |
| `process_node` | (unknown) | 24 | 25 | Near limit |

---

## References

### Code Files
- `src/tunacode/core/agents/resume/sanitize.py` - Cleanup loop with redundant canonicalization
- `src/tunacode/core/agents/resume/prune.py` - Token counting in pruning
- `src/tunacode/core/agents/main.py` - Main agent loop and logging
- `src/tunacode/core/agents/agent_components/tool_executor.py` - Tool execution with sequential batching
- `src/tunacode/core/agents/agent_components/orchestrator/tool_dispatcher.py` - Tool dispatch (McCabe 29)
- `src/tunacode/tools/parsing/command_parser.py` - JSON parsing with retry
- `src/tunacode/tools/parsing/tool_parser.py` - Fallback tool parser
- `src/tunacode/core/agents/agent_components/response_state.py` - Threading locks in asyncio
- `src/tunacode/core/types/tool_registry.py` - Tool registry with O(N) removal
- `src/tunacode/core/agents/resume/summary.py` - Unused compaction code
- `src/tunacode/core/types/state_structures.py` - Unused total_tokens field

### GitHub Permalinks
- [sanitize.py](https://github.com/alchemiststudiosDOTai/tunacode/blob/4fa9eb15d09e1263d0bd443f4ac7469b52d9904b/src/tunacode/core/agents/resume/sanitize.py) - Cleanup loop
- [prune.py](https://github.com/alchemiststudiosDOTai/tunacode/blob/4fa9eb15d09e1263d0bd443f4ac7469b52d9904b/src/tunacode/core/agents/resume/prune.py) - Pruning logic
- [tool_executor.py](https://github.com/alchemiststudiosDOTai/tunacode/blob/4fa9eb15d09e1263d0bd443f4ac7469b52d9904b/src/tunacode/core/agents/agent_components/tool_executor.py) - Tool execution
- [tool_dispatcher.py](https://github.com/alchemiststudiosDOTai/tunacode/blob/4fa9eb15d09e1263d0bd443f4ac7469b52d9904b/src/tunacode/core/agents/agent_components/orchestrator/tool_dispatcher.py) - Tool dispatch

### Delta Cards
- `.claude/delta/2026-02-02-remove-toolbuffer-shim.md` - ToolBuffer removal context
- `.claude/delta/2026-02-02-remove-heuristic-token-counting.md` - Token counting refactor

### Previous Research
- `memory-bank/research/2026-02-01_complexity-hotspots.md` - Complexity analysis
- `memory-bank/research/bottleneck-analysis.md` - Previous bottleneck analysis
