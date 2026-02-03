---
title: Decoupled Tools from Pydantic-AI
link: decoupled-tools-pydantic-ai
type: delta
path: memory-bank/delta/2026-01-28_decoupled-tools.md
depth: 1
seams: [A] architecture
ontological_relations:
  - relates_to: [[src/tunacode/tools]]
  - affects: [[src/tunacode/exceptions.py]]
  - fixes: [[direct-framework-dependency-in-tools]]
tags:
  - refactor
  - decoupling
  - architecture
created_at: 2026-01-28
---

# Summary
Removed direct dependencies on `pydantic-ai` from the `src/tunacode/tools/` directory (specifically `ModelRetry`). Tools now raise a domain-specific `ToolRetryError`, which is intercepted by the tool decorator and translated to `ModelRetry` for the framework.

# Context
Previously, 8 tool files imported `ModelRetry` from `pydantic_ai.exceptions`. This coupled the business logic of the tools directly to the specific agent framework we are using.

# Changes
- **Added `ToolRetryError`** to `src/tunacode/exceptions.py`.
- **Updated `src/tunacode/tools/decorators.py`** to catch `ToolRetryError` and raise `ModelRetry`. This file now acts as the Adapter.
- **Updated 7 tools** to import and raise `ToolRetryError` instead of `ModelRetry`:
  - `bash.py`
  - `glob.py`
  - `grep.py`
  - `list_dir.py`
  - `update_file.py`
  - `web_fetch.py`
  - `write_file.py`

# Behavioral Impact
No external behavioral change. The LLM still receives the retry signal exactly as before. The change is purely internal structural improvement.
