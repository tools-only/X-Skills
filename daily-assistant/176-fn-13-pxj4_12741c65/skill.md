# fn-13-pxj.4 Add flowctl task reset command

## Description
Add `flowctl task reset <task_id> [--cascade]` to reset task status to todo.

### Implementation

**In flowctl.py:**

1. Add `cmd_task_reset()`:
   ```python
   def cmd_task_reset(args):
       task = load_task(args.task_id)
       epic = load_epic(task["epic"])
       
       # Validations
       if epic["status"] == "done":
           error_exit("Cannot reset task in closed epic", args.json)
       if task["status"] == "in_progress":
           error_exit("Cannot reset in_progress task", args.json)
       if task["status"] == "todo":
           # Already pending, no-op success
           json_output({"success": True, "reset": [], "message": "Already todo"}) if args.json else print("Already todo")
           return
       
       # Reset status
       task["status"] = "todo"
       task["updated_at"] = now_iso()
       
       # Clear optional fields if present
       task.pop("blocked_reason", None)
       task.pop("completed_at", None)
       
       # Clear claim fields (CRITICAL for multi-user)
       task.pop("assignee", None)
       task.pop("claimed_at", None)
       task.pop("claim_note", None)
       
       # Clear evidence from JSON (if present)
       task.pop("evidence", None)
       
       save_task(task)
       
       # Clear evidence CONTENTS from task spec markdown (keep heading!)
       clear_task_evidence(args.task_id)
       
       reset_ids = [args.task_id]
       
       if args.cascade:
           # Find and reset dependent tasks (same epic only)
           dependents = find_dependents(args.task_id, same_epic=True)
           for dep_id in dependents:
               dep_task = load_task(dep_id)
               if dep_task["status"] == "in_progress":
                   continue  # Skip in_progress dependents
               if dep_task["status"] == "todo":
                   continue  # Already todo
               dep_task["status"] = "todo"
               dep_task["updated_at"] = now_iso()
               dep_task.pop("blocked_reason", None)
               dep_task.pop("completed_at", None)
               dep_task.pop("assignee", None)
               dep_task.pop("claimed_at", None)
               dep_task.pop("claim_note", None)
               dep_task.pop("evidence", None)
               save_task(dep_task)
               clear_task_evidence(dep_id)
               reset_ids.append(dep_id)
       
       # Output
       json_output({"success": True, "reset": reset_ids}) if args.json else print(f"Reset: {', '.join(reset_ids)}")
   ```

2. Add `clear_task_evidence()` helper:
   ```python
   def clear_task_evidence(task_id: str) -> None:
       """Clear ## Evidence section contents (but KEEP the heading)."""
       spec_path = get_task_spec_path(task_id)
       if not spec_path.exists():
           return
       content = spec_path.read_text()
       
       # Replace contents under ## Evidence with empty template, keeping heading
       # Pattern: ## Evidence\n<content until next ##> -> ## Evidence\n- Commits:\n- Tests:\n- PRs:\n
       # Use regex to find ## Evidence section and replace its content
       import re
       pattern = r'(## Evidence\n).*?(?=\n## |\Z)'
       replacement = r'\1- Commits:\n- Tests:\n- PRs:\n'
       new_content = re.sub(pattern, replacement, content, flags=re.DOTALL)
       
       if new_content != content:
           spec_path.write_text(new_content)
   ```

3. Add `find_dependents()` helper:
   ```python
   def find_dependents(task_id: str, same_epic: bool = False) -> list[str]:
       """Find tasks that depend on task_id (recursive)."""
       # Scan all tasks, find those with task_id in depends_on
       # If same_epic, filter to same epic
       # Recursively find their dependents too
   ```

4. Register parser:
   ```python
   p_task_reset = task_sub.add_parser("reset", help="Reset task to todo")
   p_task_reset.add_argument("task_id", help="Task ID")
   p_task_reset.add_argument("--cascade", action="store_true", help="Also reset dependents")
   p_task_reset.add_argument("--json", action="store_true")
   p_task_reset.set_defaults(func=cmd_task_reset)
   ```

### Key Behavior Notes
- `clear_task_evidence()` keeps the `## Evidence` heading (required for future `cmd_done` patches)
- Replaces evidence content with default empty template (`- Commits:\n- Tests:\n- PRs:`)
- Also clears `task["evidence"]` from JSON

### Key Files
- `plugins/flow-next/scripts/flowctl.py`
## Acceptance
- [ ] `flowctl task reset <id>` changes done/blocked → todo
- [ ] Clears optional fields: blocked_reason, completed_at (if present)
- [ ] Clears claim fields: assignee, claimed_at, claim_note
- [ ] Clears ## Evidence section from task spec markdown
- [ ] Sets updated_at to current time
- [ ] `--cascade` also resets dependent tasks (same epic only)
- [ ] Cascade skips in_progress dependents
- [ ] Errors on in_progress task
- [ ] Errors on task in closed epic
- [ ] No-op (success) if task already todo
- [ ] `--json` outputs valid JSON with reset task IDs
## Done summary
- Added `clear_task_evidence()` helper to reset spec Evidence section
- Added `find_dependents()` for recursive dependency traversal
- Added `cmd_task_reset()` with --cascade support
- Clears: status→todo, evidence, claim fields, blocked_reason

Why: Enable retry/re-run of completed tasks without manual JSON editing

Verification: Manual tests for reset, already-todo, in_progress error; CI 31/31
## Evidence
- Commits: d8fe7a45132450286a9e6777f3bc98e8874c9b1f
- Tests: plugins/flow-next/scripts/ci_test.sh
- PRs: