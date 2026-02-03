---
title: "Sanitize Loop Canonicalization Bottleneck - Plan"
phase: Plan
date: "2026-02-02"
owner: "agent"
parent_research: "memory-bank/research/2026-02-02_sanitize_loop_canonicalization_bottleneck.md"
git_commit_at_plan: "4fa9eb15"
tags: [plan, performance, canonicalization, sanitize]
---

## Goal

Reduce message canonicalization overhead in `run_cleanup_loop()` by caching canonical messages per iteration: 1N conversions when no mutations occur, and 2N conversions in iterations that mutate messages (recompute after mutation).

**Non-goals:**
- Changing the cleanup logic or iteration cap
- Modifying CanonicalMessage structure
- Adding new cleanup functions

## Scope & Assumptions

### In Scope
- `adapter.py:309-318` - `find_dangling_tool_calls()` double conversion
- `sanitize.py:469-491` - Cache canonical messages at iteration start
- Pass pre-computed canonical list to cleanup functions

### Out of Scope
- Performance of individual cleanup operations
- Changes to `to_canonical()` implementation
- CanonicalMessage class modifications

### Assumptions
- CanonicalMessage instances are immutable within an iteration
- Messages may be deleted during cleanup, invalidating cache indices; cache must be recomputed after mutation
- Current API signatures can be extended with optional canonical cache parameter

### Constraints
- Must maintain backward compatibility with existing function signatures
- No new dependencies
- Must pass existing tests

## Deliverables (DoD)

1. **Modified `find_dangling_tool_calls()`** - Single conversion per message (2N to N)
2. **Cache injection in `run_cleanup_loop()`** - Convert once at iteration start, recompute after mutation (1N no-mutation / 2N with mutation)
3. **Updated cleanup functions** - Accept optional canonical cache parameter
4. **1 test** - Verify conversion count reduction in no-mutation iterations

## Readiness (DoR)

- [x] Research document completed with line-level analysis
- [x] Code verified to match research (no drift)
- [x] Git branch ready (`bottlenecks`)
- [x] tk CLI available

## Milestones

- **M1:** Fix double conversion in `find_dangling_tool_calls` (adapter.py)
- **M2:** Add cache injection at iteration start (sanitize.py)
- **M3:** Test and verify reduction

## Work Breakdown (Tasks)

### Task 1: Fix Double Conversion in find_dangling_tool_calls
**Owner:** agent
**Dependencies:** None
**Milestone:** M1

**What:**
Modify `adapter.py:309-318` to convert each message once and extract both call_ids and return_ids from the same canonical instance.

**Before:**
```python
for msg in messages:
    call_ids.update(get_tool_call_ids(msg))     # to_canonical(msg)
    return_ids.update(get_tool_return_ids(msg)) # to_canonical(msg) AGAIN
```

**After:**
```python
for msg in messages:
    canonical = msg if isinstance(msg, CanonicalMessage) else to_canonical(msg)
    call_ids.update(canonical.get_tool_call_ids())
    return_ids.update(canonical.get_tool_return_ids())
```

**Files:**
- `src/tunacode/utils/messaging/adapter.py:309-318`

**Acceptance:**
- [ ] Each message converted at most once per call
- [ ] Function returns same results as before
- [ ] Existing tests pass

---

### Task 2: Add Canonical Cache at Iteration Start
**Owner:** agent
**Dependencies:** Task 1
**Milestone:** M2

**What:**
At the start of each cleanup iteration in `run_cleanup_loop()`, convert messages once and pass the canonical list to `find_dangling_tool_call_ids()`. If any cleanup mutates `messages`, recompute the canonical list before subsequent cleanup steps.

**Files:**
- `src/tunacode/core/agents/resume/sanitize.py:469-491`

**Acceptance:**
- [ ] `to_canonical_list()` called once at iteration start
- [ ] Canonical list passed to `find_dangling_tool_call_ids()`
- [ ] Canonical cache recomputed after any mutation before later cleanup passes

---

### Task 3: Update Cleanup Functions for Cache Passthrough
**Owner:** agent
**Dependencies:** Task 2
**Milestone:** M2

**What:**
Modify `remove_empty_responses()` and `remove_consecutive_requests()` to accept an optional canonical cache parameter, skipping their internal `_canonicalize_messages()` call when cache is provided.

**Files:**
- `src/tunacode/core/agents/resume/sanitize.py:326-355` (remove_empty_responses)
- `src/tunacode/core/agents/resume/sanitize.py:358-404` (remove_consecutive_requests)

**Acceptance:**
- [ ] Functions accept optional `canonical_cache` parameter
- [ ] When cache provided, internal conversion skipped
- [ ] Existing tests pass

---

### Task 4: Write Verification Test
**Owner:** agent
**Dependencies:** Task 3
**Milestone:** M3

**What:**
Add one test that verifies conversion count reduction in a no-mutation iteration. Use a mock or counter to track `to_canonical()` calls.

**Files:**
- `tests/unit/test_sanitize_canonicalization.py` (new)

**Acceptance:**
- [ ] Test verifies N conversions per iteration when no mutations occur (1N, not 4N)
- [ ] Test passes

## Risks & Mitigations

| Risk | Impact | Likelihood | Mitigation | Trigger |
|------|--------|------------|------------|---------|
| Cache recompute missed after mutation | High | Medium | Recompute canonical cache after any deletion before later cleanup passes | Index errors or mismatched removals |
| Backward compatibility break | Medium | Low | Keep original signatures, add optional params | Import errors in other modules |
| Test flakiness with mock counters | Low | Low | Use deterministic test data | Random failures in CI |

## Test Strategy

- **1 new test:** `test_canonicalization_count_reduction` - Verify conversion count is 1N (not 4N) in a no-mutation iteration using a mock counter on `to_canonical()`.

## References

- Research: `memory-bank/research/2026-02-02_sanitize_loop_canonicalization_bottleneck.md`
- `adapter.py:309-318` - Double conversion bottleneck
- `sanitize.py:469-491` - Cleanup loop

## Tickets Created

| Ticket ID | Title | Priority | Status |
|-----------|-------|----------|--------|
| bot-555f | Fix double conversion in find_dangling_tool_calls | P1 | open |
| bot-30b4 | Add canonical cache at iteration start in run_cleanup_loop | P1 | open |
| bot-11c6 | Update cleanup functions for cache passthrough | P2 | open |
| bot-a9de | Write verification test for conversion count reduction | P2 | open |

## Dependencies

```
bot-555f (ready)
    └── bot-30b4 (blocked by bot-555f)
        └── bot-11c6 (blocked by bot-30b4)
            └── bot-a9de (blocked by bot-11c6)
```
