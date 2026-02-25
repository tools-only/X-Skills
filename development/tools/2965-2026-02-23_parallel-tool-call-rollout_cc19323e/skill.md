---
type: delta
date: 2026-02-23
status: done
tags: [tinyagent, parallel-tool-calls, concurrency-cap, tui]
---

# Parallel Tool-call Rollout via tiny-agent-os 1.2.5

## Scope completed

- `tun-d506` — bumped `tiny-agent-os` to `>=1.2.5` and refreshed `uv.lock`.
- `tun-q44p` — enforced a hard cap of 3 concurrent tool executions by wrapping all runtime `AgentTool.execute` handlers with a shared semaphore in `agent_config.py`.
- `tun-5nny` — hardened `RequestOrchestrator` tool lifecycle handling:
  - strict `tool_execution_start` args normalization (`dict` or JSON object string only)
  - fail-loud `TypeError` on malformed payloads
  - explicit duration semantics for parallel batches (suppress per-tool duration when batch execution is detected)
- `tun-6qpn` — updated TUI status behavior for parallel starts/completions:
  - status bar tracks multiple in-flight tools coherently
  - callbacks now emit completion before last-action updates
  - shell status path aligned with the same completion semantics

## Repo truth updates

- README migration note now reflects parallel batches with max 3 in-flight calls.
- Changelog `Unreleased` section documents dependency cutover, runtime cap, orchestrator hardening, and UI behavior update.

## Validation

- `uv run ruff check --fix .`
- `uv run ruff check .`
- `uv run pytest tests/unit/core/test_tool_concurrency_limit.py tests/unit/core/test_request_orchestrator_parallel_tools.py tests/unit/ui/test_status_bar_layout.py tests/system/cli/test_repl_support.py tests/unit/ui/test_tool_panel_css_flow.py`
- `uv run pytest`
- `uv run python -c "import importlib.metadata as m; print(m.version('tiny-agent-os'))"` -> `1.2.5`
