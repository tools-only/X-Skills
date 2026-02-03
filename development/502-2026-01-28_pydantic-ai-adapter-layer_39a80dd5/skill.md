---
title: "Decouple Tools from Pydantic-AI"
phase: Completed
date: "2026-01-28"
owner: "architecture"
tags: [architecture, decoupling, pydantic-ai, refactor]
---

## Goal
Remove the direct dependency on `pydantic-ai` from `src/tunacode/tools/`.
Tools should raise domain-specific exceptions, not framework-specific exceptions.

## Strategy
1.  **Define Protocol:** Create `ToolRetryError` in `src/tunacode/exceptions.py`.
2.  **Adapter Layer:** Update `src/tunacode/tools/decorators.py` to catch `ToolRetryError` and re-raise `pydantic_ai.ModelRetry`.
3.  **Refactor Tools:** Update all tools to raise `ToolRetryError` instead of `ModelRetry`.

## Steps (Executed 2026-01-28)

### Step 1: Define Exception
- Added `ToolRetryError` to `src/tunacode/exceptions.py`.

### Step 2: Update Decorator
- Modified `src/tunacode/tools/decorators.py` to catch `ToolRetryError` and re-raise `ModelRetry`.
- This acts as the translation layer between our internal tool protocol and the Pydantic-AI framework.

### Step 3: Migrate Tools
- Replaced `ModelRetry` with `ToolRetryError` in:
    - `tools/bash.py`
    - `tools/glob.py`
    - `tools/grep.py`
    - `tools/list_dir.py`
    - `tools/update_file.py`
    - `tools/web_fetch.py`
    - `tools/write_file.py`

## Verification
`grep -r "pydantic_ai" src/tunacode/tools/` returns only `decorators.py`.
Dependencies are now correctly directed: `tools` -> `exceptions` <- `decorator` -> `pydantic-ai`.
