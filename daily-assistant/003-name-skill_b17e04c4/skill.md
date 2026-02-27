---
name: orchestrator-discipline
description: Orchestrator context window discipline enforcement. Prevents the orchestrator from reading source files it will not edit, running diagnostic commands that waste context, and rationalizing delegation bypasses. Use when setting up orchestrator guardrails, reviewing delegation discipline, or diagnosing context window waste in multi-agent workflows. Activates PreToolUse hooks that surface decision points before source file reads and diagnostic command execution.
user-invocable: true
---
# Orchestrator Discipline

Enforce delegation discipline for the orchestrator role in multi-agent Claude Code workflows. The orchestrator's context window is a shared resource across the entire session — agents get fresh context per task, the orchestrator does not.

## What This Plugin Provides

### 1. PreToolUse Hooks (Structural Enforcement)

Two hooks fire automatically on every tool call:

**Source File Read Warning** — fires on `Read` and `Grep` targeting source/config/test files (`.py`, `.toml`, `.yaml`, `.json`, etc.). Injects a decision-point reminder asking: "Will you Edit/Write this file this turn?"

**Diagnostic Command Gate** — fires on `Bash` calls matching diagnostic commands (`ty check`, `ruff check`, `mypy`, `pytest`, `eslint`, `cargo check`, etc.). Reminds the orchestrator to delegate the command to an Explore agent instead.

Both hooks are **non-blocking** — they inject `additionalContext` to surface the decision, not prevent it. Legitimate reads (reading a file you are about to edit) proceed normally.

### 2. Rules (Behavioral Constraints)

The `rules/CLAUDE.md` file is loaded into every session and provides:

- Read permission/prohibition lists with a falsifiable test
- Delegation constraint definitions (no exemption categories)
- Investigation escalation anti-pattern documentation
- Diagnostic command delegation patterns
- Epistemic identity scoping for orchestrator role

### 3. Reference Material

See [Investigation Escalation Anti-Pattern](./references/investigation-escalation.md) for the detailed pattern analysis, root cause diagnosis, and correct workflow alternatives.

## When to Activate

This skill auto-loads via the plugin's hooks and rules. Manual activation is useful when:

- Reviewing whether the orchestrator is following delegation discipline
- Diagnosing context window waste in a session
- Training new orchestrator configurations on delegation patterns

## Correct Orchestrator Workflow

```mermaid
flowchart TD
    Start([Task arrives]) --> Q1{Does orchestrator need<br>current codebase state?}
    Q1 -->|Yes| Explore["Delegate to Explore agent<br>with diagnostic command"]
    Q1 -->|No| Direct[Delegate implementation<br>to specialist agent]
    Explore --> Summary[Receive summary<br>from agent]
    Summary --> Q2{Scope changed?}
    Q2 -->|Yes| User[Present to user<br>for routing decision]
    Q2 -->|No| Direct
    User --> Direct
    Direct --> PostCheck["Spot-check agent output<br>(read deliverable only)"]

    Q1 -.->|WRONG| Self["Read source files yourself<br>Run diagnostics yourself"]
    Self -.->|Escalates to| Bypass["'This is simple enough<br>to do myself'"]
```

## Hook Behavior Reference

### Source File Read Warning

**Triggers on**: `Read` or `Grep` where target path matches:

- Extensions: `.py`, `.toml`, `.yaml`, `.yml`, `.js`, `.ts`, `.jsx`, `.tsx`, `.json`, `.cfg`, `.ini`, `.env`, `.sh`, `.bash`, `.go`, `.rs`, `.rb`, `.java`, `.c`, `.cpp`, `.h`, `.hpp`
- Test paths: directories named `test/`, `tests/`, `spec/`, `__tests__/`, or files matching `test_*.py`

**Does NOT trigger on**: `.md`, `.txt`, plan files, backlog items, CLAUDE.md, skill definitions

### Diagnostic Command Gate

**Triggers on**: `Bash` where command matches:

- Python: `ty check`, `ruff check`, `mypy`, `pyright`, `basedpyright`, `pylint`, `pytest`
- JavaScript/TypeScript: `eslint`, `tsc --noEmit`
- Rust: `cargo check`, `cargo clippy`
- Go: `go vet`
- Meta: `pre-commit run`, `prek run`

**Does NOT trigger on**: `git status`, `ls`, `wc`, `uv run` (without diagnostic subcommand), or any non-diagnostic bash command
