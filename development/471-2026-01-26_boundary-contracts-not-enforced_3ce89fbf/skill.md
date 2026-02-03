---
title: Boundary contracts are implicit across callbacks and tools
link: boundary-contracts-not-enforced
type: delta
path: src/tunacode/types/callbacks.py
depth: 0
seams: [A, M]
ontological_relations:
  - relates_to: [[architecture]]
  - affects: [[ui-core-boundary]]
  - affects: [[core-tools-boundary]]
tags:
  - architecture
  - boundaries
  - callbacks
  - protocols
created_at: 2026-01-26T03:09:30Z
updated_at: 2026-01-26T03:09:30Z
uuid: 278ae3c6-3c42-4c1c-94ff-b935c148aa21
ticket: t-69f0
---

## Summary

Boundary contracts between UI, core, and tools are partially implicit: several callbacks and tool entry points are typed as `Any` or use broad protocols, so the contract is not enforced by typing or tests. This makes boundary drift harder to detect and allows accidental coupling across layers.

## Context

Recent refactors established clean import direction and layer ordering, but contract surfaces (callbacks, tool factories, and state accessors) remain loosely typed in a few places.

## Root Cause

Contracts were defined informally via usage instead of explicit callback type aliases and minimal protocols. We missed it because there is no boundary-contract audit or type gate for these callback signatures; prevention is to encode contracts in `types/` and add a focused type/check step.

## Changes

- Pending: define explicit callback type aliases for tool results, plan approval, and streaming.
- Pending: narrow `StateManagerProtocol` into per-tool protocols (todo, plan approval, authorization).
- Pending: update tool factories and UI hooks to accept only the minimal protocols and callbacks.

## Behavioral Impact

**What users notice:**
- No immediate behavior change, but boundary leaks can manifest as subtle UI/core coupling and harder-to-debug callback misuse.

**What didn't change:**
- Runtime behavior and tool execution flow remain unchanged.

## Related Cards

- [[tool-registry-todo-alignment]]
