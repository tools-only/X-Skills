---
name: implement-feature
description: Execute a SAM task plan (plan/tasks-*.md) by looping ready tasks, delegating each to its specified agent, and relying on hooks to update task timestamps/status. Use when a task file exists and you need to run the implementation loop that picks up ready tasks and delegates them to agents.
argument-hint: <task-file-path or feature-slug>
user-invocable: true
hooks:
  SubagentStop:
  - hooks:
    - type: command
      command: python3 "${CLAUDE_PLUGIN_ROOT}/skills/implementation-manager/scripts/task_status_hook.py"
version: 1.0.0
last_updated: '2026-01-27'
python_compatibility: 3.11+
---
# Implement Feature (SAM Workflow Execution)

This workflow continues from `add-new-feature`. It executes tasks from a SAM task file until complete (or blocked).

<feature_input>
$ARGUMENTS
</feature_input>

---

## Resolve Task File

Rules:

- If `$ARGUMENTS` ends with `.md`, treat it as the task file path.
- Otherwise, treat it as a feature slug (or partial slug) and resolve via `implementation_manager.py`.

Example resolution:

```bash
uv run "${CLAUDE_PLUGIN_ROOT}/skills/implementation-manager/scripts/implementation_manager.py" \
  status "${CLAUDE_PROJECT_DIR}" "$ARGUMENTS"
```

---

## Progress Loop

1. Query status:

```bash
uv run "${CLAUDE_PLUGIN_ROOT}/skills/implementation-manager/scripts/implementation_manager.py" \
  status "${CLAUDE_PROJECT_DIR}" "$ARGUMENTS"
```

2. If tasks remain, query ready tasks:

```bash
uv run "${CLAUDE_PLUGIN_ROOT}/skills/implementation-manager/scripts/implementation_manager.py" \
  ready-tasks "${CLAUDE_PROJECT_DIR}" "$ARGUMENTS"
```

3. For each ready task:

- Route to the agent named in the task's `**Agent**` field.
- Launch the agent with a prompt that invokes `start-task`:

```text
Skill(skill="python3-development:start-task", args="{task_file_path} --task {task_id}")
```

4. Repeat until no tasks remain ready.

---

## Completion Gate

When all tasks show `COMPLETE`, invoke:

```text
Skill(skill="python3-development:complete-implementation", args="{task_file_path}")
```
