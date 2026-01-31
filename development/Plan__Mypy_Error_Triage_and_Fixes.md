---
name: Plan: Mypy Error Triage and Fixes
source: https://raw.githubusercontent.com/alchemiststudiosDOTai/tunacode/master/PLAN.md
original_path: PLAN.md
source_repo: alchemiststudiosDOTai/tunacode
category: development
subcategory: coding
tags: ['development']
collected_at: 2026-01-31T18:34:05.962039
file_hash: e5c04865865de2cbda3ec8de7002e31966bd31cd4443a124c6a0b8f539ed18ce
---

# Plan: Mypy Error Triage and Fixes

## Objective
Reduce mypy errors incrementally by investigating root causes first, then fixing in small, safe batches.

## Current State
- Last full mypy run (before Batch 3/4 updates) reported 19 errors in 7 files.
- **Action needed:** re-run full mypy to refresh the error list after Batches 3–5.

## Investigation Phase
- ✅ Verified initial error list and classified by category/module.

## Fix Plan (Incremental Batches)

### Batch 1: Low-risk typing hygiene ✅
**Scope:** missing return/argument/variable annotations.
**Files updated:**
- `src/tunacode/tools/utils/ripgrep.py`
- `src/tunacode/tools/parsing/retry.py`
- `src/tunacode/utils/system/gitignore.py`
- `src/tunacode/configuration/settings.py`
- `src/tunacode/configuration/paths.py`
- `src/tunacode/tools/grep.py`
- `src/tunacode/tools/glob.py`
- `src/tunacode/tools/bash.py`
- `src/tunacode/core/state.py`
**Verification:** mypy for Batch 1 files passed.

### Batch 2: Optional handling & assignment mismatches ✅
**Scope:** optional cache handling + ripgrep mismatch.
**Changes:**
- Fixed None-attribute access in `src/tunacode/tools/ignore.py` with early return.
- Ripgrep mismatch not reproducible in current code.
**Verification:** mypy on `ignore.py` and `ripgrep.py` passed.

### Batch 3: Protocol/interface mismatches ✅
**Scope:** `StateManagerProtocol` alignment + orchestration callback signature.
**Changes:**
- Switched agent components to `StateManagerProtocol`.
- Expanded `SessionStateProtocol` for fields used by agents.
- Fixed orchestrator tool_result_callback invocation signature.
**Files updated:**
- `src/tunacode/core/agents/agent_components/agent_config.py`
- `src/tunacode/core/agents/agent_components/streaming.py`
- `src/tunacode/core/agents/agent_components/orchestrator/orchestrator.py`
- `src/tunacode/core/agents/agent_components/orchestrator/tool_dispatcher.py`
- `src/tunacode/core/agents/agent_components/orchestrator/message_recorder.py`
- `src/tunacode/core/agents/agent_components/orchestrator/usage_tracker.py`
- `src/tunacode/core/types/state.py`
**Verification:** mypy for touched files passed.

### Batch 4: UI contract mismatches ✅
**Scope:** UI protocols and callback signatures.
**Changes:**
- Aligned `ShellRunnerHost.notify` with Textual notify signature.
- Updated `AppForCallbacks.post_message` return type and `status_bar` typing.
- Fixed `_on_setup_complete`, TUI callback typing, and `_request_worker` return type.
**Files updated:**
- `src/tunacode/ui/shell_runner.py`
- `src/tunacode/ui/repl_support.py`
- `src/tunacode/ui/app.py`
**Verification:** mypy for UI files passed.

### Batch 5: Coroutine/task API correctness ✅
**Scope:** task handling in grep.
**Verification:** mypy for `src/tunacode/tools/grep.py` passed (no code changes needed).

## Quality Gates (last run)
- `ruff check --fix /home/tuna/tunacode` ✅
- `uv run mypy /home/tuna/tunacode/src/tunacode` ❌ (before Batch 3/4 updates)
- `uv run pytest` ✅ (510 passed)
- `grimp` dependency check ✅ (0 violations)

## Next Steps
1. Re-run full mypy to refresh remaining errors.
2. Create new batch plan for any remaining files.

## Definition of Done (per batch)
- Mypy errors reduced for the batch scope only.
- No new type errors introduced.
- Minimal diffs with explicit typing and early returns.

## Notes
- Use small, focused diffs.
- Avoid cross-layer dependency violations.
- Update docs only if behavior changes.
