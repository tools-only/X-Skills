# @File Passing Root-Fix Plan (No Cutovers)

## Hypothesis (Current Best)

The primary break is **our @file message pipeline**, not tinyagent/alchemy typing.

Evidence in current code:

- `src/tunacode/ui/widgets/file_autocomplete.py`
  - completion inserts values that still include `@` (e.g. `@pyproject.toml`), not absolute paths.
- `src/tunacode/ui/app.py`
  - `on_editor_submit_requested` enqueues raw `message.text` directly (`request_queue.put(message.text)`).
- `src/tunacode/ui/repl_support.py`
  - no active `@` normalization helper in HEAD.

So we currently pass unresolved mention tokens into the agent loop.

---

## Direction Lock

- Keep tinyagent/alchemy typed tool-call flow unchanged.
- No cutover shims.
- Fix @mention ingestion at the UI boundary.
- Add evidence logging to confirm payload shape after @ fix.

---

## Scope

### In-scope

1. `src/tunacode/ui/repl_support.py`
2. `src/tunacode/ui/app.py`
3. `tests/system/cli/test_repl_support.py`
4. `tests/unit/ui` (new focused submit-path test)
5. (optional evidence) `src/tunacode/tools/decorators.py` warning log only

### Out-of-scope (for this pass)

- Tool registry ID remapping
- adapter/sanitize ID/name fallback logic
- type alias removals

We only widen scope if evidence still shows malformed tool calls after @ fix.

-- Implementation Plan

### 1) Add deterministic @mention normalization

**File:** `src/tunacode/ui/repl_support.py`

Add `normalize_agent_message_text(text: str, *, cwd: Path | None = None) -> str`:

- Detect tokens of form `@<path>` bounded by whitespace/start.
- Strip trailing punctuation (`.,:;!?)]}`) for path resolution, then re-attach.
- Resolve relative paths against `cwd` (default `Path.cwd()`), expanduser, absolute resolve.
- Replace only when resolved path exists.
- Leave unresolved mentions untouched (fail loud later via tool result, no silent guessing).

Rule: no aliasing, no magic; strict path existence check.

---

### 2) Wire normalization into submit path

**File:** `src/tunacode/ui/app.py`

In `on_editor_submit_req quested(...)`:

- Keep command handling on raw submitted text.
- Normalize text once before queuing request.
- Queue normalized text.
- Keep rendered user message as raw typed text.

This preserves UX while fixing execution payload.

---

### 3) Add evidence logging for chain verification

**Files:**

- `src/tunacode/ui/app.py` (submit boundary)
- optionally `src/tunacode/tools/decorators.py` (tool execute boundary)

Logging (warning-level only while investigating):

- raw submitted text (redacted/truncated)
- normalized text (redacted/truncated)
- whether replacement occurred
- at tool execute: args keys/type before bind

Goal: prove whether malformed tool args persist after @ normalization.

---

### 4) Tests

#### A. `tests/system/cli/test_repl_support.py`

Add cases:

1. `@relative/path` resolves to absolute existing path.
2. punctuation preserved (`@file.py,` -> `/abs/file.py,`).
3. non-existent mention unchanged.
4. message without `@` unchanged.

#### B. New unit UI submit-path test

**Suggested file:** `tests/unit/ui/test_app_editor_submit_mentions.py`

Validate:

- raw text is what user sees.
- queued text is normalized text.
- `/command` handling still uses raw text and short-circuits queueing.

---

## Verification Commands

- `uv run pytest tests/system/cli/test_repl_support.py -q`
- `uv run pytest tests/unit/ui/test_app_editor_submit_mentions.py -q`
- `uv run ruff check src/tunacode/ui/repl_support.py src/tunacode/ui/app.py tests/system/cli/test_repl_support.py tests/unit/ui/test_app_editor_submit_mentions.py`

---

## Acceptance Criteria

1. `@file` mentions are converted to absolute existing paths before agent request execution.
2. User-visible chat text remains unchanged from user input.
3. Command parsing behavior is unchanged.
4. Evidence logs can confirm whether malformed tool args still occur post-fix.
5. No tinyagent/alchemy contract changes introduced.

---

## Escalation Rule

If malformed tool payloads still appear after this fix (confirmed by logs), then open Phase 2 for runtime hardening in registry/adapter/sanitize. Not before.
