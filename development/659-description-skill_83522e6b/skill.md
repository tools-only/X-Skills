---
description: Orchestrator delegation workflows for linting. Guides orchestrators on when and how to delegate to linting-root-cause-resolver and post-linting-architecture-reviewer agents. Use when orchestrating linting tasks, delegating quality checks, or reading linting resolution reports.
---

# Holistic Linting: Orchestrator Delegation

This skill provides orchestrator-specific workflows for delegating linting and resolution tasks to specialized agents.

## When to Use This Skill

Use this skill when you are an **orchestrator** (interactive Claude Code CLI) after completing implementation work. This skill guides delegation patterns for linting, quality checks, and architectural review.

**Do NOT use this skill** if you are a sub-agent executing a delegated task. Sub-agents should follow the linter-specific resolution workflows in the holistic-linting-resolver skill.

## Agent Delegation (Orchestrator Only)

### Complete Linting Workflow

**CRITICAL PRINCIPLE**: Orchestrators delegate work to agents. Orchestrators do NOT run formatting commands, linting commands, or quality checks themselves. The agent does ALL work (formatting, linting, resolution). The orchestrator only delegates tasks and reads reports to determine if more work is needed.

**WHY THIS MATTERS**:

- Pre-gathering linting data wastes orchestrator context window
- Running linters duplicates agent work (agent will run them again)
- Violates separation of concerns: "Orchestrators route context, agents do work"
- Creates context rot with linting output that becomes stale
- Prevents agents from gathering their own fresh context

The orchestrator MUST follow this delegation-first workflow:

**Step 1: Delegate to linting-root-cause-resolver immediately**

Delegate linting resolution WITHOUT running any linting commands first:

```text
Task(
  agent="holistic-linting:linting-root-cause-resolver",
  prompt="Format, lint, and resolve any issues in <file_path>"
)
```

**What NOT to do before delegating**:

- ❌ Do NOT run `ruff format` before delegating
- ❌ Do NOT run `ruff check` before delegating
- ❌ Do NOT run `mypy` before delegating
- ❌ Do NOT gather linting output for the agent
- ❌ Do NOT read error messages to provide to the agent

**What TO do**:

- ✅ Delegate immediately with just the file path
- ✅ Let agent gather its own linting data
- ✅ Trust agent to run formatters and linters itself
- ✅ Wait for agent to complete and produce reports

**Reason**: The agent follows systematic root-cause analysis workflows. It autonomously:

- Discovers project linters by scanning configuration files
- Runs formatters on modified files (ruff format, prettier, etc.)
- Executes linters to identify issues (ruff, mypy, pyright, etc.)
- Researches rule documentation
- Traces type flows and architectural context
- Implements elegant fixes following python3-development patterns
- Verifies resolution by re-running linters
- Creates resolution artifacts in `.claude/reports/` and `.claude/artifacts/`

**Multiple Files Modified**:

Launch concurrent agents (one per file) WITHOUT pre-gathering linting data:

```text
Task(agent="holistic-linting:linting-root-cause-resolver", prompt="Format, lint, and resolve any issues in src/auth.py")
Task(agent="holistic-linting:linting-root-cause-resolver", prompt="Format, lint, and resolve any issues in src/api.py")
Task(agent="holistic-linting:linting-root-cause-resolver", prompt="Format, lint, and resolve any issues in tests/test_auth.py")
```

**Reason for concurrency**: Independent file resolutions proceed in parallel, reducing total time.

**Step 2: Delegate to post-linting-architecture-reviewer**

After linting agent completes, delegate architectural review:

```text
Task(
  agent="post-linting-architecture-reviewer",
  prompt="Review linting resolution for <file_path>"
)
```

**What the reviewer does**:

- Loads resolution artifacts from `.claude/reports/` and `.claude/artifacts/`
- Verifies resolution quality (root cause addressed, no symptom suppression)
- Validates architectural implications (design principles, type safety, code organization)
- Identifies systemic improvements applicable across codebase
- Generates architectural review report

**Step 3: Read reviewer report**

The orchestrator reads the review report to determine if additional work is needed:

```bash
ls -la .claude/reports/architectural-review-*.md
```

Read the most recent review report:

```claude
Read(".claude/reports/architectural-review-[timestamp].md")
```

**Orchestrator's role**: Read reports and decide next steps. Do NOT run linting commands to verify agent's work.

**Step 4: If issues found, delegate back to linting agent**

If architectural review identifies problems with resolution:

```text
Task(
  agent="holistic-linting:linting-root-cause-resolver",
  prompt="Address issues found in architectural review: .claude/reports/architectural-review-[timestamp].md

Issues identified:
- [Summary of finding 1]
- [Summary of finding 2]

Review report contains detailed context and proposed solutions."
)
```

**Step 5: Repeat review if needed**

After re-resolution, delegate to reviewer again:

```text
Task(
  agent="post-linting-architecture-reviewer",
  prompt="Review updated linting resolution for <file_path>"
)
```

Continue workflow until architectural review reports clean results.

### Workflow Summary

```text
[Implementation complete]
  → [Step 1: Delegate to linting-root-cause-resolver] (agent formats, lints, resolves)
  → [Step 2: Delegate to post-linting-architecture-reviewer]
  → [Step 3: Orchestrator reads review report]
  → [Step 4: If issues found, delegate back to resolver with review path]
  → [Step 5: Repeat review until clean]
  → [Task complete ✓]
```

**Key Principle**: Orchestrator delegates immediately and reads reports. Agent does ALL actionable work (formatting, linting, resolution). Orchestrator does NOT run commands or gather linting data.

### Common Anti-Patterns to Avoid

**❌ WRONG** - Orchestrator pre-gathering linting data:

```text
# Don't do this:
Bash("ruff check src/auth.py")
# Read the output...
# Then delegate with the output
Task(agent="holistic-linting:linting-root-cause-resolver", prompt="Fix these errors: [pasted errors]")
```

**✅ CORRECT** - Orchestrator delegates immediately:

```text
# Do this instead:
Task(agent="holistic-linting:linting-root-cause-resolver", prompt="Format, lint, and resolve any issues in src/auth.py")
```

**❌ WRONG** - Orchestrator running formatters:

```text
# Don't do this:
Bash("ruff format src/auth.py src/api.py")
# Then delegate linting
```

**✅ CORRECT** - Agent handles both formatting and linting:

```text
# Do this instead:
Task(agent="holistic-linting:linting-root-cause-resolver", prompt="Format, lint, and resolve any issues in src/auth.py")
```

**❌ WRONG** - Orchestrator verifying agent's work by running linters:

```text
# Don't do this:
# Agent completes
Bash("ruff check src/auth.py")  # Verifying agent's work
```

**✅ CORRECT** - Trust agent's verification, read reports instead:

```text
# Do this instead:
Read(".claude/reports/linting-resolution-[timestamp].md")
# Report shows agent already verified with linter output
```

## Related Skills

- [holistic-linting](../holistic-linting/SKILL.md) - Core linting skill with linter detection and resource documentation
- [holistic-linting-resolver](../holistic-linting-resolver/SKILL.md) - Linter-specific resolution workflows for sub-agents
