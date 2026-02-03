---
title: Dangling Tool Calls on User Abort
link: dangling-tool-calls-abort
type: delta
path: src/tunacode/core/agents/main.py
depth: 0
seams: [M, S]
ontological_relations:
  - relates_to: [[conversation-turns]]
  - affects: [[message-history]]
  - fixes: [[api-error-on-next-request]]
tags:
  - bug
  - abort-handling
  - message-invariants
  - tool-calls
created_at: 2026-01-17T06:00:00Z
updated_at: 2026-01-17T06:00:00Z
uuid: d7c8e9f0-1a2b-3c4d-5e6f-7a8b9c0d1e2f
---

# Dangling Tool Calls on User Abort

## Summary

When a user aborted mid-tool-call (Ctrl+C or tool denial), the conversation history was left with a `ModelResponse` containing tool calls but no corresponding `ToolReturn` messages. The next API request failed because the provider expected tool returns for pending calls.

## Context

- File: `src/tunacode/core/agents/main.py`
- PR: #246
- Surfaces: On any abort during tool execution, next request fails

## Root Cause

**Missing invariant enforcement on exception path.**

The message flow has an implicit invariant:

```
INVARIANT: Every ModelResponse with tool_calls MUST be followed by
           matching ToolReturn(s) before the next ModelRequest
```

The happy path naturally maintains this - tool executes, result recorded. But `UserAbortError` breaks out of the loop BEFORE recording tool returns, violating the invariant.

```
┌─────────────────────────────────────────────────────────────┐
│ Happy Path:                                                 │
│   ModelResponse(tool_calls=[A, B])                          │
│   → execute tools                                           │
│   → ToolReturn(A), ToolReturn(B)                            │
│   → ModelRequest (next turn)  ✓                             │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│ Abort Path (BEFORE fix):                                    │
│   ModelResponse(tool_calls=[A, B])                          │
│   → execute tool A                                          │
│   → UserAbortError raised                                   │
│   → messages still has ModelResponse with [A, B]            │
│   → NO ToolReturn recorded                                  │
│   → ModelRequest (next turn)  ✗ API ERROR                   │
└─────────────────────────────────────────────────────────────┘
```

## Changes

Added `_remove_dangling_tool_calls()` in the `except UserAbortError` handler:

```python
except UserAbortError:
    cleanup_applied = _remove_dangling_tool_calls(
        self.state_manager.session.messages,
        self.state_manager.session.tool_call_args_by_id,
    )
    if cleanup_applied:
        self.state_manager.session.update_token_count()
    raise
```

The function walks backward through messages, removing any trailing `ModelResponse` that has unanswered tool calls, and clears cached args.

## Behavioral Impact

- User aborts mid-tool → dangling messages cleaned up
- Next request works normally
- Token count recalculated after cleanup

## Prevention: What We Should Have Done

### 1. Document Message Invariants

The conversation turns doc now exists, but we missed documenting this invariant:

```markdown
## Message Invariants

1. **Tool Call Pairing**: Every `ModelResponse` with tool_calls MUST
   be followed by matching `ToolReturn(s)` before any new `ModelRequest`.

2. **Exception Safety**: Any exception path that exits the agent loop
   must restore message history to a valid state.
```

Add this to `docs/codebase-map/architecture/conversation-turns.md`.

### 2. Add Exception Path Tests

No test existed for abort scenarios:

```python
async def test_abort_mid_tool_call_cleans_up_messages():
    """Verify abort during tool execution doesn't corrupt message history."""
    # Arrange: Start a request that will trigger a tool call
    # Act: Raise UserAbortError during tool execution
    # Assert: messages ends with a valid state (no dangling tool calls)
    # Assert: Next request succeeds
```

### 3. State Machine Should Track Message Validity

The `ResponseState` tracks iteration progress but not message-level invariants. Consider:

```python
def validate_message_invariants(messages: list[Any]) -> bool:
    """Return False if messages violate tool-call pairing invariant."""
    # Walk through messages
    # Track pending tool calls
    # Verify each is answered before next ModelRequest
```

This could be called:
- After each agent iteration (debug mode)
- Before persisting session
- On session load

### 4. Design Pattern: Exception-Safe State Updates

The broader lesson: **state mutations during a loop must be exception-safe**.

Options:
1. **Rollback on exception** (what we did) - detect and fix invalid state
2. **Transactional updates** - only commit state after successful completion
3. **Copy-on-write** - work on a copy, swap on success

For message history, option 2 would be:
```python
# Don't append to session.messages during the loop
# Collect in a local buffer
# Merge to session.messages only on successful completion
```

This is partially what `agent_run.all_messages()` does - but the exception path didn't use it.

## Related Cards

- [[conversation-turns]] - Full turn flow documentation
- [[message-history-persistence]] - How messages are persisted
- [[response-state]] - State machine for agent progress

## Lesson

**Exception paths are first-class citizens.** When designing a stateful loop:

1. List every exception that can exit the loop
2. For each: what state is left behind?
3. For each: is that state valid for the next operation?
4. Add tests for exception scenarios

The bug was obvious in hindsight - we just never traced what happens when `UserAbortError` interrupts the tool execution flow.
