---
title: Parallel tool-call rollout plan for TunaCode
type: metadata
link: parallel-tool-calls-rollout-plan
path: .claude/metadata/research/2026-02-23-parallel-tool-calls-rollout-plan.md
depth: 2
seams: [A]
ontological_relations:
  - affects: [[core-agents]]
  - affects: [[ui]]
  - affects: [[tools]]
  - relates_to: [[tinyagent]]
tags:
  - tinyagent
  - parallel-tool-calls
  - rollout-plan
  - tui
created_at: 2026-02-23T14:54:29-06:00
updated_at: 2026-02-23T14:54:29-06:00
uuid: ac4317c5-b747-4efe-bd6f-039f4bdc3833
---

# Parallel tool-call rollout plan for TunaCode

## Objective

Adopt tinyagent parallel tool-call execution in TunaCode with no compatibility shim, while preserving correctness, UI clarity, and fail-loud behavior.

## Code-state snapshot (from current repo)

- Dependency floor is still `tiny-agent-os>=1.2.4` (`pyproject.toml:31`).
- Lockfile currently resolves `tiny-agent-os` to `1.2.4` (`uv.lock` package stanza).
- README still states sequential execution (`README.md:40`).
- TunaCode runtime already consumes tinyagent tool lifecycle events:
  - `tool_execution_start` handler (`src/tunacode/core/agents/main.py:396`)
  - `tool_execution_end` handler (`src/tunacode/core/agents/main.py:417`)
- UI callbacks are start/end based and status-bar driven:
  - `build_tool_start_callback()` (`src/tunacode/ui/repl_support.py:220`)
  - `build_tool_result_callback()` (`src/tunacode/ui/repl_support.py:189`)
  - status methods (`src/tunacode/ui/widgets/status_bar.py:52`, `:55`)
- Vendored tinyagent code is already on 1.2.5 semantics:
  - parallel gather in `execute_tool_calls()` (`tinyAgent/tinyagent/agent_tool_execution.py:186`)
  - all starts emitted before ends (`tinyAgent/tests/test_parallel_tool_execution.py`)
  - post-batch steering semantics documented in `tinyAgent/CHANGELOG.md`.

## Implementation plan

### Phase 1 — Dependency cutover (no shim)

1. Bump TunaCode dependency floor to `tiny-agent-os>=1.2.5`.
2. Refresh lockfile with `uv lock` so resolved version is >=1.2.5.
3. Verify runtime import path remains unchanged (TunaCode already uses tinyagent public API).

**Acceptance**
- `pyproject.toml` and `uv.lock` both reflect 1.2.5-compatible resolution.
- `uv run python -c "import importlib.metadata as m; print(m.version('tiny-agent-os'))"` prints >=1.2.5 in project env.

### Phase 2 — Runtime hardening for parallel batches

1. **Normalize tool args at event ingress** in `RequestOrchestrator` start-handler.
   - Replace cast-only behavior with explicit type checks.
   - If args are invalid (non-dict JSON object), fail loud with a clear error.
2. **Validate duration semantics** under batch completion ordering.
   - Current per-tool duration is measured start→end event in TunaCode.
   - In parallel mode, all end events are emitted after gather completion, so durations may converge to batch latency.
   - Decide one policy and apply consistently:
     - either keep current metric but label it as batch-relative,
     - or suppress per-tool duration when multiple tools are active in a batch.
3. Keep steering logic aligned to tinyagent contract (post-batch application).

**Acceptance**
- No callback crashes from malformed `args` payloads.
- Tool duration presentation is intentional and documented in code comments/tests.

### Phase 3 — UI behavior alignment for concurrent starts

1. Verify status bar behavior when multiple start events arrive before end events.
2. Keep existing minimal behavior (`running: <tool>`, then `last: <tool>`) for first cut, or upgrade to explicit concurrent state (count/list) in same PR if low risk.
3. Ensure tool panels remain classed correctly (`running/completed/failed`) through `tool_panel_smart()`.

**Acceptance**
- No UI regressions in tool panel rendering.
- Status bar stays coherent during multi-tool turns.

### Phase 4 — Tests and gates

Add focused tests around TunaCode event handling (not tinyagent internals):

1. `RequestOrchestrator` handles parallel-style order:
   - start(A), start(B), end(A), end(B)
   - registry and callbacks remain correct per tool.
2. `RequestOrchestrator` arg-normalization path:
   - invalid args payload triggers deterministic failure.
3. Optional UI test for status bar transitions under batched starts.

Run gates:
- `uv run pytest`
- `uv run ruff check --fix .`
- `uv run ruff .`

### Phase 5 — Repo truth updates

1. Update README statement that still claims sequential-only tool execution.
2. Add changelog entry describing the cutover and any intentional UI semantics.

## Risks to actively watch

1. **Concurrent mutating tools** (e.g., two file writes/edits in one assistant turn) can create race-like behavior in shared file/cache state.
2. **Metric trust** if per-tool durations look precise but are batch-relative.
3. **User visibility** during concurrent tool starts if status bar only shows one tool name.

## Suggested rollout order

1. Dependency bump + lockfile
2. Runtime hardening + tests
3. UI polish (if needed)
4. README/changelog updates
