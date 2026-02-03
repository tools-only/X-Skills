---
title: Restore tunacode.cli.main Import Compatibility
link: cli-main-compatibility
type: delta
path: src/tunacode/cli/main.py
depth: 0
seams: [E, M]
ontological_relations:
  - relates_to: [[cli-entrypoint]]
  - affects: [[python-package]]
  - fixes: [[cli-import-compatibility]]
tags:
  - cli
  - compatibility
  - release
  - lint
created_at: 2026-02-02T12:45:00-06:00
updated_at: 2026-02-02T12:45:00-06:00
uuid: c7977f36-a39f-4e9b-b972-d93c1a7a359b
---

# Restore tunacode.cli.main Import Compatibility

## Summary

Release validation exposed that `tunacode.cli.main` no longer existed after the CLI entrypoint moved to `tunacode.ui.main`, breaking import checks and any downstream scripts using the old path. We added a small compatibility shim that re-exports the Typer app and aligned linted `isinstance` checks with PEP 604 unions to keep release tooling green.

## Context

- File: `src/tunacode/cli/main.py`
- Symptom: `ModuleNotFoundError` on `from tunacode.cli.main import app` during release testing.

## Root Cause

The CLI entrypoint was consolidated under `tunacode.ui.main` without leaving a compatibility package for the previous import path, so older references failed at import time.

## Changes

- Added a `tunacode.cli` package with a `main.py` shim that re-exports `app` from `tunacode.ui.main`.
- Re-exported `app` in `tunacode.cli.__init__` for convenience.
- Swapped tuple-based `isinstance` checks to PEP 604 union syntax to satisfy ruff pre-flight checks.

## Behavioral Impact

- `tunacode.cli.main` imports now resolve without errors.
- CLI behavior remains unchanged because the entrypoint still delegates to `tunacode.ui.main`.

## Related Cards

- [[cli-entrypoint]]
