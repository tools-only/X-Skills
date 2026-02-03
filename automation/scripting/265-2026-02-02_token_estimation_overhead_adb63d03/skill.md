# Research: Token Estimation Overhead in Pruning

**Date:** 2026-02-02
**Phase:** Research
**Scope:** `src/tunacode/core/agents/resume/prune.py:161-238`

## Structure

### Token Estimation Pipeline

```
RequestOrchestrator.run()
  └─> _prepare_message_history()           [main.py:278]
      └─> _log_pruned_tool_outputs()       [main.py:310]
          └─> prune_old_tool_outputs()     [prune.py:161]
              ├─> estimate_part_tokens()   [prune.py:100]    <-- BOTTLENECK
              │   └─> estimate_tokens()    [token_counter.py:10]
              └─> prune_part_content()     [prune.py:121]
                  └─> estimate_tokens()    [token_counter.py:10]  (redundant)
```

### Directory Layout

```
src/tunacode/
├── utils/messaging/
│   ├── token_counter.py    # Core heuristic: len(text) // 4
│   └── __init__.py         # Exports estimate_tokens()
└── core/agents/resume/
    ├── prune.py            # Backward-scan algorithm (161-238)
    ├── summary.py          # Compaction trigger (uses estimate_tokens)
    ├── sanitize.py         # Cleanup loop (separate bottleneck)
    └── filter.py           # Calls prune_old_tool_outputs()
```

## Key Files

### Core Algorithm

- `prune.py:161-238` → `prune_old_tool_outputs()` - main entry point
- `prune.py:100-118` → `estimate_part_tokens()` - per-part token estimation
- `prune.py:121-158` → `prune_part_content()` - mutation and reclaim calculation

### Token Estimation

- `token_counter.py:7` → `CHARS_PER_TOKEN = 4` constant
- `token_counter.py:10-14` → `estimate_tokens()` function: `len(text) // 4`

### Callers

- `main.py:310-312` → `_log_pruned_tool_outputs()` triggers pruning
- `filter.py:74` → `prepare_history()` triggers pruning

## Algorithm Phases (prune.py:161-238)

### Phase 1: Collection (Lines 188-203)
```python
# O(M × P × L) - scans ALL messages backward
for msg_idx in range(len(messages) - 1, -1, -1):
    for part_idx in range(len(parts) - 1, -1, -1):
        if is_tool_return_part(part):
            tokens = estimate_part_tokens(part, model_name)  # O(L)
            tool_parts.append((msg_idx, part_idx, part, tokens))
```

### Phase 2: Boundary (Lines 210-218)
```python
# O(T) - finds 40k token protection boundary
for i, (_, _, _, tokens) in enumerate(tool_parts):
    accumulated_tokens += tokens
    if accumulated_tokens > protect_tokens:  # 40,000
        prune_start_index = i
        break
```

### Phase 3: Savings (Lines 224-230)
```python
# O(T) - validates 20k minimum threshold
parts_to_prune = tool_parts[prune_start_index:]
total_prunable_tokens = sum(tokens for _, _, _, tokens in parts_to_prune)
if total_prunable_tokens < minimum_threshold:  # 20,000
    return (messages, 0)
```

### Phase 4: Mutation (Lines 232-238)
```python
# O(T × L) - REDUNDANT token estimation
for _, _, part, _ in parts_to_prune:
    reclaimed = prune_part_content(part, model_name)  # re-estimates tokens!
```

## Complexity Analysis

### Time Complexity: O(M × P × L)

| Variable | Description | Typical Range |
|----------|-------------|---------------|
| M | Number of messages | 50-500 |
| P | Parts per message | 1-5 |
| L | Content length (chars) | 100-50,000 |

**Worst case example:**
- 200 messages × 3 parts × 10,000 chars = 6,000,000 character scans

### Space Complexity: O(M × P)

- `tool_parts` list stores up to M × P tuples
- No deep copies (references only)

## Patterns Found

### Pattern 1: Redundant Token Calculation

**Location:** Phase 1 calculates tokens (line 201), Phase 4 recalculates (lines 144, 149)

```python
# Phase 1 - already calculated
tokens = estimate_part_tokens(part, model_name)
tool_parts.append((msg_idx, part_idx, part, tokens))  # stored here!

# Phase 4 - recalculated from scratch
for _, _, part, _ in parts_to_prune:  # token count in tuple IGNORED
    reclaimed = prune_part_content(part, model_name)
    # prune_part_content calls estimate_tokens() twice more (lines 144, 149)
```

### Pattern 2: No Memoization

- Same content estimated multiple times
- No caching by content hash or part ID
- `PRUNE_PLACEHOLDER` tokens recalculated per part (line 149)

### Pattern 3: Unused Parameter

- `model_name` passed through all functions
- Never actually used in estimation (heuristic is model-agnostic)

## Dependencies

### Import Chain

```
estimate_tokens() defined in:
  └─ utils/messaging/token_counter.py:10

Re-exported through:
  └─ utils/messaging/__init__.py:16

Imported by:
  ├─ core/agents/resume/prune.py:11
  └─ core/agents/resume/summary.py:19
```

### Caller Chain

```
main.py:278 _prepare_message_history()
  └─> main.py:310 _log_pruned_tool_outputs()
      └─> prune.py:161 prune_old_tool_outputs()

filter.py:74 prepare_history()
  └─> prune.py:161 prune_old_tool_outputs()
```

## Symbol Index

### Constants (prune.py)

| Symbol | Line | Value |
|--------|------|-------|
| `PRUNE_PROTECT_TOKENS` | 14 | 40,000 |
| `PRUNE_MINIMUM_THRESHOLD` | 15 | 20,000 |
| `PRUNE_MIN_USER_TURNS` | 17 | 2 |
| `PRUNE_PLACEHOLDER` | 18 | "[Old tool result content cleared]" |
| `PART_KIND_TOOL_RETURN` | 37 | "tool-return" |
| `PART_KIND_USER_PROMPT` | 38 | "user-prompt" |

### Constants (token_counter.py)

| Symbol | Line | Value |
|--------|------|-------|
| `CHARS_PER_TOKEN` | 7 | 4 |

### Functions (prune.py)

| Function | Line | Signature |
|----------|------|-----------|
| `get_prune_thresholds` | 27 | `() -> tuple[int, int]` |
| `is_tool_return_part` | 41 | `(part: Any) -> bool` |
| `is_user_prompt_part` | 59 | `(part: Any) -> bool` |
| `count_user_turns` | 73 | `(messages: list[Any]) -> int` |
| `estimate_part_tokens` | 100 | `(part: Any, model_name: str) -> int` |
| `prune_part_content` | 121 | `(part: Any, model_name: str) -> int` |
| `prune_old_tool_outputs` | 161 | `(messages: list[Any], model_name: str) -> tuple[list[Any], int]` |

## Data Structures

### tool_parts List (prune.py:190)

```python
tool_parts: list[tuple[int, int, Any, int]]
#                      │    │    │    └─ token_count (calculated)
#                      │    │    └────── part object (reference)
#                      │    └─────────── part_index in message.parts
#                      └──────────────── message_index in messages list
```

**Order:** Most recent first (backward scan appends newest to end, but iteration order is newest-first)

## Bottleneck Analysis

### Primary Bottleneck: Phase 1 Collection

**Why it's expensive:**
1. Iterates ALL messages regardless of early exit conditions
2. Calls `estimate_part_tokens()` for EVERY tool return part
3. Each call does `len(content)` on potentially large strings

**Measured complexity:**
- 200 messages × 3 tool parts × 5000 char average = 3M character operations
- Called on every agent request

### Secondary Bottleneck: Phase 4 Redundancy

**Why it's wasteful:**
1. Token counts already in `tool_parts` tuple (index 3)
2. `prune_part_content()` re-estimates original tokens (line 144)
3. `prune_part_content()` re-estimates placeholder tokens (line 149)
4. Placeholder is constant - should be pre-calculated once

### Optimization Opportunities

1. **Use pre-calculated tokens from Phase 1** in Phase 4 instead of re-estimating
2. **Cache PRUNE_PLACEHOLDER token count** as module constant
3. **Early termination** in Phase 1 once enough tokens collected (don't scan entire history)
4. **Batch token estimation** - concatenate contents, single length call
5. **Incremental tracking** - store token counts on parts, update only when content changes
