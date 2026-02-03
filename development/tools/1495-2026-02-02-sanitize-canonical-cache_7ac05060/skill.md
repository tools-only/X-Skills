---
title: Cache Canonical Messages in Cleanup Loop
link: sanitize-canonical-cache
type: delta
path: src/tunacode/core/agents/resume/sanitize.py
depth: 0
seams: [M, D]
ontological_relations:
  - relates_to: [[message-history]]
  - affects: [[sanitize-loop]]
  - fixes: [[canonicalization-bottleneck]]
tags:
  - performance
  - sanitize
  - canonicalization
  - cleanup-loop
created_at: 2026-02-03T03:20:03Z
updated_at: 2026-02-03T03:20:03Z
uuid: ac2bdaf8-9aae-4027-8841-170a46503d96
---

# Cache Canonical Messages in Cleanup Loop

## Summary

The sanitize cleanup loop re-canonicalized messages multiple times per iteration. We now cache canonical messages per iteration, recomputing only after mutations, and remove the double conversion inside `find_dangling_tool_calls()`.

## Context

- Files: `src/tunacode/core/agents/resume/sanitize.py`, `src/tunacode/utils/messaging/adapter.py`
- Symptom: redundant `to_canonical()` calls during cleanup iterations.

## Root Cause

Canonicalization was happening in each cleanup helper (and twice per message in `find_dangling_tool_calls()`), so a single iteration performed multiple full conversions of the same message list.

## Changes

- Converted each message once inside `find_dangling_tool_calls()` by reusing a single canonical instance.
- Cached canonical messages at the start of each cleanup iteration and recomputed only after mutations.
- Added canonical cache passthroughs to `remove_empty_responses()` and `remove_consecutive_requests()`.
- Added a unit test to verify 1N conversions in a no-mutation iteration.

## Behavioral Impact

- Cleanup iterations now run with fewer canonicalization passes (1N without mutations, 2N when deletions occur).
- Cleanup behavior and tool-call pruning logic remain unchanged.

## Related Cards

- [[message-history]]
- [[sanitize-loop]]
