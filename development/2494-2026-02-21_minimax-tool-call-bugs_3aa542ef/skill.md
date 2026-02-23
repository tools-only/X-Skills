---
type: delta
date: 2026-02-21
status: open
tags: [minimax, tool-call, bug]
---

# MiniMax Tool Call Bugs -- Handoff

## Root Cause

MiniMax M2.5 returns tool calls with empty `id` and sometimes empty `name` fields. Our Python code rejects empty strings in three places. The tinyagent Rust layer (alchemy-llm) handles IDs correctly via `ToolCallId` with `Into` trait -- this is NOT a tinyagent bug.

## Three rejection sites (all ours)

1. **`src/tunacode/core/types/tool_registry.py:44`** -- `register()` raises `ValueError("tool_call_id is required")` on empty string
2. **`src/tunacode/utils/messaging/adapter.py:155`** -- `_coerce_tool_call_item()` raises `TypeError` on empty `id`
3. **`src/tunacode/core/agents/resume/sanitize.py:147`** -- raises `TypeError` on empty `id`

Fix: generate fallback UUID at each site instead of raising. Registry fix requires adding `tool_id_map: dict[str, str]` to `_TinyAgentStreamState` (`main.py:75`) so the end handler resolves fallback IDs.

## `read_file` argument mismatch

MiniMax ignores the JSON schema parameter name `filepath`. Unknown what key it actually sends -- debug logging attempt failed because logger level was DEBUG (filtered out). Use `logging.warning` to capture actual `args` dict in `decorators.py` `execute()` before `sig.bind()`.

Once you know the actual key, add it to the alias list in `_normalize_tool_args`.

System prompt examples should use JSON object format (`read_file({"filepath": "..."})`) not function-call syntax, since that's what the model actually emits.

## `@file` resolution

Not yet implemented. Add `normalize_agent_message_text` to `repl_support.py` to resolve `@path` tokens to absolute paths. Wire into `app.py` `on_editor_submit_requested`. Note: `@` triggers editor autocomplete, so tmux send-keys won't work for testing -- test manually or via unit test.

## Dependency

Bump `tiny-agent-os>=1.2.3` in `pyproject.toml` (separate from these bugs, but needed for other fixes in that release).

## Stash

`git stash list` has a previous attempt. Partially usable code for `decorators.py` normalization and `repl_support.py` @file resolution. The rest was wrong.
