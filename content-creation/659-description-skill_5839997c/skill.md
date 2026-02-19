---
description: Start or complete a specific task inside a SAM task file. Updates task status to IN PROGRESS with Started timestamp, writes active-task context for hooks, and supports --complete to mark tasks complete. Use when an agent needs to pick up a SAM task, set it in progress, implement acceptance criteria, and signal completion.
argument-hint: <task-file-path> [--task <task-id>] [--complete <task-id>]
user-invocable: true
hooks:
  PostToolUse:
  - matcher: Write|Edit|Bash
    hooks:
    - type: command
      command: python3 "${CLAUDE_PLUGIN_ROOT}/skills/implementation-manager/scripts/task_status_hook.py"
version: 1.0.0
last_updated: '2026-01-27'
python_compatibility: 3.11+
---

# Start Task (SAM Task Execution Helper)

You are implementing a specific task from a SAM task file.

<task_input>
$ARGUMENTS
</task_input>

---

## Parse Arguments

- `task_file_path` (required): path to a `plan/tasks-*.md` file
- `--task <id>` (optional): Task ID to start (defaults to first ready task)
- `--complete <id>` (optional): Task ID to mark COMPLETE

---

## Detect Task Format

Read the task file. Determine the format:

- **YAML frontmatter**: File starts with `---` and contains YAML fields like `task:`, `status:`, `title:`
- **Individual task file**: A single file with YAML frontmatter representing ONE task
- **Inline markdown**: File contains `## Task N:` headers with `**Status**:` bold fields

For individual task files (YAML frontmatter at top), the task IS the entire file.
For monolithic files (multiple tasks), find the specific task section.

---

## If `--complete <task-id>` Provided

1. Read the task file.
2. Update the selected task:

   **If YAML frontmatter format:**
   - Edit `status:` field to `complete`
   - Add `completed: {ISO timestamp}` field

   **If inline markdown format:**
   - Change `**Status**` to `✅ COMPLETE`
   - Add/update `**Completed**: {ISO timestamp}`

3. Output: `Task {ID} marked as complete`

---

## Starting a Task

1. Read the task file and the linked architecture spec.
2. Select the task:
   - If `--task` provided, use that ID
   - Else pick the first task where status is `not-started` (YAML) or `NOT STARTED` (markdown) and all dependencies are resolved
3. Update the task status:

   **If YAML frontmatter format:**
   - Edit the `status:` field in frontmatter to `in-progress`
   - Add `started: {ISO timestamp}` field to frontmatter

   **If inline markdown format:**
   - Set `**Status**: 🔄 IN PROGRESS`
   - Add `**Started**: {ISO timestamp}`

4. Write the active-task context file (required for hook-driven updates):

```bash
mkdir -p .claude/context
printf '%s' '{"task_file_path": "{task_file_path}", "task_id": "{task_id}"}' > ".claude/context/active-task-${CLAUDE_SESSION_ID}.json"
```

5. Implement against the task acceptance criteria and run its verification steps.
