---
title: UI token estimation routed through core messaging facade
type: delta
link: ui-messaging-layer-facade
path: src/tunacode/ui/app.py
depth: 0
seams: [A, M]
ontological_relations:
  - relates_to: [[dependency-layers]]
  - affects: [[ui]]
  - affects: [[core-ui-api]]
  - fixes: [[ui-imports-utils-layer-violation]]
tags:
  - architecture
  - ui
  - dependencies
  - messaging
created_at: 2026-02-05T15:56:00-06:00
updated_at: 2026-02-05T15:56:00-06:00
uuid: 9f8848d1-f855-4567-99a8-fa37a8f104de
---

# UI token estimation routed through core messaging facade

## Summary

Fixed a dependency-layer violation where the UI imported `tunacode.utils.messaging` directly by routing message token estimation through the `tunacode.core.ui_api.messaging` facade.

## Context

`tests/test_dependency_layers.py::test_no_new_layer_violations` flagged:

- `tunacode.ui.app -> tunacode.utils.messaging` (ui -> utils)

## Root Cause

The UI layer is only allowed to depend on `core` (per the dependency layer rules), but `TextualReplApp._update_resource_bar()` imported `estimate_messages_tokens` from `utils`.

## Changes

- `src/tunacode/core/ui_api/messaging.py`
  - Re-exported `estimate_messages_tokens()` from `tunacode.utils.messaging`.
- `src/tunacode/ui/app.py`
  - Switched import to `from tunacode.core.ui_api.messaging import estimate_messages_tokens`.
  - Reordered imports to satisfy first-party import-order tests.

## Behavioral Impact

No user-visible behavior change; token estimation remains the same. The UI now respects the enforced dependency direction.

## Related Cards

- [[dependency-layers]]
