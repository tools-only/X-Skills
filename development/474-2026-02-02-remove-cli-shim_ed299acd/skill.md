---
title: Remove tunacode.cli Compatibility Shim
link: remove-cli-shim
type: delta
path: src/tunacode/cli
depth: 0
seams: [E, M]
ontological_relations:
  - relates_to: [[cli-entrypoint]]
  - affects: [[python-package]]
  - fixes: [[cli-import-compatibility]]
tags:
  - cli
  - compatibility
  - removal
created_at: 2026-02-02T12:43:36-06:00
updated_at: 2026-02-02T12:43:36-06:00
uuid: 7e644a16-2214-4da2-acc7-7c1fcd834519
---

# Remove tunacode.cli Compatibility Shim

## Summary

Removed the temporary `tunacode.cli` package shim so the only supported CLI entrypoint is `tunacode.ui.main`, aligning code and docs on a single import path.

## Context

- Removed compatibility package introduced to satisfy a release-time import check.
- Updated module documentation to reflect the canonical package path.

## Root Cause

The shim was a short-term workaround and no longer reflects the desired single-source CLI entrypoint.

## Changes

- Deleted `src/tunacode/cli/` package.
- Updated command parser module docstring to the correct package path.
- Added changelog note under Unreleased about the removal.

## Behavioral Impact

- Importing `tunacode.cli.main` will now fail; use `tunacode.ui.main` instead.
- CLI behavior remains unchanged because the runtime entrypoint already uses `tunacode.ui.main`.

## Related Cards

- [[cli-entrypoint]]
