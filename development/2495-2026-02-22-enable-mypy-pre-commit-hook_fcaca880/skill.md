---
title: Enable mypy at pre-commit stage (keep pre-push enforcement)
type: delta
link: enable-mypy-pre-commit-hook
path: .pre-commit-config.yaml
depth: 0
seams: [S]
ontological_relations:
  - affects: [[tooling]]
  - affects: [[typing-contracts]]
  - affects: [[developer-workflow]]
tags:
  - mypy
  - pre-commit
  - quality-gate
created_at: 2026-02-22T15:17:04-06:00
updated_at: 2026-02-22T15:17:04-06:00
uuid: 6a5c9232-30f7-4aa9-bf94-ecb0197d6263
---

# Enable mypy at pre-commit stage (keep pre-push enforcement)

## Summary

Updated `.pre-commit-config.yaml` so the `mypy` hook runs during both:

- `pre-commit`
- `pre-push`

This catches type errors earlier at commit time while retaining the existing push-time gate.

## Verification

- `uv run pre-commit run mypy --files src/tunacode/ui/widgets/editor.py` ✅
