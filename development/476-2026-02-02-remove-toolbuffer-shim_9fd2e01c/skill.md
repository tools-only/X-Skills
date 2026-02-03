---
title: Remove ToolBuffer shim and strengthen tool batching prompt
link: remove-toolbuffer-shim-prompt-batching
type: delta
path: src/tunacode/core/agents
depth: 1
seams: [M]
ontological_relations:
  - relates_to: [[core-agents]]
  - affects: [[system-prompt]]
tags:
  - core
  - agents
  - tools
  - prompt
  - cleanup
created_at: 2026-02-02T14:19:29-06:00
updated_at: 2026-02-02T14:19:29-06:00
uuid: e56b5cc0-6e0f-4b52-8d97-f7c831219734
---

## Summary
Removed the unused ToolBuffer shim from the agent loop and strengthened the system prompt to demand batched read-only tool calls during discovery, aligning the runtime with actual parallel tool execution capabilities.

## Context
ToolBuffer had no call sites and was only flushed at request completion, so it never influenced execution. Prompt guidance needed explicit batching language to reduce sequential read loops.

## Changes
- Removed ToolBuffer file, exports, and call wiring from the agent loop and orchestrator signature.
- Updated parallel execution guidance in the system prompt with explicit batching requirements and an example.

## Behavioral Impact
- No runtime behavior change from ToolBuffer removal (buffer was unused).
- Model guidance now more forcefully requests batched tool calls after discovery.

## Related Cards
- [[2026-02-02-tool-name-normalization]]
