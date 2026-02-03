---
title: Pruning Token Estimation Optimization
link: pruning-token-estimation-optimization
type: delta
path: src/tunacode/core/agents/resume/prune.py
depth: 2
seams: [M]
ontological_relations:
  - relates_to: [[prune_old_tool_outputs]]
  - relates_to: [[token_estimation]]
  - affects: [[context-window]]
tags:
  - performance
  - pruning
  - tokens
created_at: 2026-02-02T21:05:00-06:00
updated_at: 2026-02-02T21:05:00-06:00
uuid: 0172237b-d7e0-40ca-b130-e8dcc5bd72d2
---

# Summary

Reduced redundant token estimation in tool-output pruning while preserving pruning behavior.

# Context

The pruning loop in `prune_old_tool_outputs()` was estimating token counts twice for pruned parts and recalculating placeholder tokens for every mutation.

# Changes

- Added `PRUNE_PLACEHOLDER_TOKENS` to cache placeholder token counts.
- Introduced `get_part_content_text()` to normalize part content before estimation/pruning.
- Passed precomputed token counts into `prune_part_content()` to avoid re-estimation.
- Stopped token estimation after `protect + minimum` is satisfied while still scanning for parts to prune.

# Behavioral Impact

Pruning decisions and output remain unchanged; token estimation work is reduced and placeholder token counting is cached.

# Related Cards

- [[token-pruning-efficiency-plan]]
