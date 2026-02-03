---
title: Normalize Tool Names Before Dispatch
link: tool-name-normalization
type: delta
path: src/tunacode/core/agents/agent_components/orchestrator/tool_dispatcher.py
depth: 0
seams: [M, D]
ontological_relations:
  - relates_to: [[tool-dispatcher]]
  - affects: [[tool-calls]]
  - fixes: [[tool-name-whitespace]]
tags:
  - bug
  - tool-dispatch
  - normalization
  - retry
created_at: 2026-02-02T02:17:43Z
updated_at: 2026-02-02T02:17:43Z
uuid: 2b3b6716-d365-492f-b14f-55e69c9cd4fe
---

# Normalize Tool Names Before Dispatch

## Summary

Tool calls with leading/trailing whitespace in `tool_name` failed to dispatch (e.g., " glob"), resulting in unknown-tool errors instead of retries. We now normalize tool names before registration and dispatch to avoid whitespace-induced failures.

## Context

- File: `src/tunacode/core/agents/agent_components/orchestrator/tool_dispatcher.py`
- Symptom: tool names like " glob" were treated as unknown tools.
- Gap: tests covered suspicious characters but not whitespace-only tool names.

## Root Cause

Tool names were used verbatim from model output, and the suspicious-name guard did not treat whitespace as invalid. This allowed malformed names to reach dispatch without correction.

## Changes

- Added `_normalize_tool_name()` to trim whitespace and fall back to `UNKNOWN_TOOL_NAME` for empty results.
- Normalized tool names during structured and fallback tool-call registration.
- Updated tool-call parts in place only when normalization changes the name.

## Behavioral Impact

- Tool calls with incidental whitespace now dispatch correctly.
- Empty/whitespace-only tool names fail fast as `unknown` and trigger model retry prompts.
- No change to tool execution retry/backoff behavior.

## Related Cards

- [[tool-dispatcher]]
- [[tool-calls]]
