# Research: Sanitize Loop Canonicalization Bottleneck

**Date:** 2026-02-02
**Phase:** Research
**Branch:** bottlenecks

## Problem Statement

`run_cleanup_loop()` performs 4N message conversions per iteration (up to 10 iterations), resulting in worst-case 40N conversions for N messages.

## Structure

### File Layout
```
src/tunacode/
├── core/agents/resume/sanitize.py     # Cleanup loop orchestrator
├── utils/messaging/adapter.py          # Message conversion + tool call extraction
└── types/canonical.py                   # CanonicalMessage definition
```

### Entry Point
- `sanitize.py:458-502` - `run_cleanup_loop(messages, tool_registry)`
- Called from: `main.py:291`

## Key Files

### sanitize.py - Cleanup Loop

| Location | Function | Purpose |
|----------|----------|---------|
| `sanitize.py:38` | `MAX_CLEANUP_ITERATIONS = 10` | Loop iteration cap |
| `sanitize.py:61-63` | `_canonicalize_messages()` | Delegates to `to_canonical_list()` |
| `sanitize.py:248-250` | `find_dangling_tool_call_ids()` | Delegates to adapter |
| `sanitize.py:326-355` | `remove_empty_responses()` | Converts all messages at line 332 |
| `sanitize.py:358-404` | `remove_consecutive_requests()` | Converts all messages at line 364 |
| `sanitize.py:458-502` | `run_cleanup_loop()` | Main loop - 3 cleanup steps per iteration |

### adapter.py - Conversion Functions

| Location | Function | Conversions |
|----------|----------|-------------|
| `adapter.py:146` | `to_canonical(message)` | 1 per call |
| `adapter.py:205-207` | `to_canonical_list(messages)` | N per call |
| `adapter.py:295-299` | `get_tool_call_ids(msg)` | 1 per call (if not already canonical) |
| `adapter.py:302-306` | `get_tool_return_ids(msg)` | 1 per call (if not already canonical) |
| `adapter.py:309-318` | `find_dangling_tool_calls(messages)` | 2N per call (see bottleneck 1) |

## Patterns Found

### Pattern 1: Double Conversion in find_dangling_tool_calls

**Location:** `adapter.py:309-318`

```python
def find_dangling_tool_calls(messages: list[Any]) -> set[str]:
    call_ids: set[str] = set()
    return_ids: set[str] = set()

    for msg in messages:                         # N iterations
        call_ids.update(get_tool_call_ids(msg))   # → to_canonical(msg)  [conversion 1]
        return_ids.update(get_tool_return_ids(msg))  # → to_canonical(msg)  [conversion 2]

    return call_ids - return_ids
```

**Observation:** Same message converted twice per loop iteration. The canonical result from `get_tool_call_ids()` is discarded before `get_tool_return_ids()` re-converts the same message.

### Pattern 2: Conversion-for-Inspection in remove_empty_responses

**Location:** `sanitize.py:326-355`

```python
def remove_empty_responses(messages: list[Any]) -> bool:
    canonical_messages = _canonicalize_messages(messages)  # Line 332: N conversions
    # ... inspection logic ...
    del messages[index]  # Line 351: mutates ORIGINAL list, not canonical
```

**Observation:** Canonical messages used only for read-only inspection. Mutations happen on original pydantic-ai list.

### Pattern 3: Identical Conversion-for-Inspection in remove_consecutive_requests

**Location:** `sanitize.py:358-404`

```python
def remove_consecutive_requests(messages: list[Any]) -> bool:
    canonical_messages = _canonicalize_messages(messages)  # Line 364: N conversions
    # ... inspection logic ...
    del messages[index]  # Line 400: mutates ORIGINAL list
```

**Observation:** Same pattern as Pattern 2. Independent conversion with no shared cache.

### Pattern 4: No Caching Between Functions

**Location:** `sanitize.py:469-491`

```python
for cleanup_iteration in range(MAX_CLEANUP_ITERATIONS):
    dangling_tool_call_ids = find_dangling_tool_call_ids(messages)  # 2N conversions
    # ...
    removed_empty_responses = remove_empty_responses(messages)       # N conversions
    # ...
    removed_consecutive_requests = remove_consecutive_requests(messages)  # N conversions
```

**Observation:** Three functions called sequentially, each converts the entire message list independently. No shared canonical cache within an iteration.

## Dependencies

### Conversion Call Chain

```
run_cleanup_loop (sanitize.py:458)
├── find_dangling_tool_call_ids (sanitize.py:248)
│   └── find_dangling_tool_calls (adapter.py:309)
│       ├── get_tool_call_ids (adapter.py:295)
│       │   └── to_canonical (adapter.py:146)  # 1 per msg
│       └── get_tool_return_ids (adapter.py:302)
│           └── to_canonical (adapter.py:146)  # 1 per msg (DUPLICATE)
├── remove_empty_responses (sanitize.py:326)
│   └── _canonicalize_messages (sanitize.py:61)
│       └── to_canonical_list (adapter.py:205)  # N conversions
└── remove_consecutive_requests (sanitize.py:358)
    └── _canonicalize_messages (sanitize.py:61)
        └── to_canonical_list (adapter.py:205)  # N conversions
```

### Import Dependencies

```
sanitize.py imports from adapter.py:
├── find_dangling_tool_calls
└── to_canonical_list (via _canonicalize_messages)

sanitize.py imports from canonical.py:
├── CanonicalMessage
├── MessageRole
└── SystemPromptPart
```

## Conversion Count Summary

### Per Iteration (N messages)

| Function | Location | Conversions |
|----------|----------|-------------|
| `find_dangling_tool_call_ids()` via `get_tool_call_ids()` | `adapter.py:315` | N |
| `find_dangling_tool_call_ids()` via `get_tool_return_ids()` | `adapter.py:316` | N |
| `remove_empty_responses()` | `sanitize.py:332` | N |
| `remove_consecutive_requests()` | `sanitize.py:364` | N |
| **Total per iteration** | | **4N** |

### Total Conversions

| Scenario | Messages | Iterations | Total Conversions |
|----------|----------|------------|-------------------|
| Best case | 100 | 1 | 400 |
| Typical | 100 | 3 | 1,200 |
| Worst case | 100 | 10 | 4,000 |
| Worst case | 500 | 10 | 20,000 |
| Worst case | 1000 | 10 | 40,000 |

## Symbol Index

### sanitize.py Exports

- `sanitize_history_for_resume()` - Main entry point (line 412)
- `run_cleanup_loop()` - Loop orchestrator (line 458)
- `find_dangling_tool_call_ids()` - ID extraction (line 248)
- `remove_dangling_tool_calls()` - Removal (line 253)
- `remove_empty_responses()` - Cleanup (line 326)
- `remove_consecutive_requests()` - Cleanup (line 358)

### adapter.py Exports (relevant)

- `to_canonical()` - Single message conversion (line 146)
- `to_canonical_list()` - Batch conversion (line 205)
- `get_tool_call_ids()` - ID extraction (line 295)
- `get_tool_return_ids()` - ID extraction (line 302)
- `find_dangling_tool_calls()` - Detection (line 309)

## Bottleneck Root Causes

1. **Double conversion in find_dangling_tool_calls**: Same message converted twice to extract `call_ids` then `return_ids`

2. **Independent conversions per function**: Each cleanup function (`find_dangling`, `remove_empty`, `remove_consecutive`) converts the full list independently

3. **No per-iteration cache**: Canonical messages created and discarded by each function; no shared cache within an iteration

4. **Read-only conversions on mutable data**: Conversions used only for inspection while mutations happen on original pydantic-ai list

## Fix Locations

| Fix | File | Lines | Impact |
|-----|------|-------|--------|
| Single conversion in find_dangling | `adapter.py` | 309-318 | -N per iteration |
| Cache canonical at iteration start | `sanitize.py` | 469-491 | -2N per iteration |
| Pass cache to cleanup functions | `sanitize.py` | 326, 364 | Eliminates redundant calls |

**After fixes:** 1N conversions per iteration (down from 4N) = 75% reduction
