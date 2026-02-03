# Plan-Sync Agent Specification

> **Implementation ready.** This spec contains all details needed for autonomous implementation.

## File Paths (all relative to `plugins/flow-next/`)

| File | Full Path |
|------|-----------|
| flowctl.py | `plugins/flow-next/scripts/flowctl.py` |
| phases.md | `plugins/flow-next/skills/flow-next-work/phases.md` |
| plan-sync.md | `plugins/flow-next/agents/plan-sync.md` (CREATE) |
| smoke_test.sh | `plugins/flow-next/scripts/smoke_test.sh` |

## Overview

A subagent that synchronizes downstream task specs when implementation drifts from the original plan. Runs after each task completes, comparing what was specified vs what was built, and updates affected downstream tasks.

**Problem solved:** Tasks planned early may reference assumptions that become stale during implementation. Variable names change, APIs evolve, approaches shift. Downstream tasks inherit outdated assumptions, causing confusion or wasted work.

**Solution:** After each task completes, plan-sync:
1. Re-anchors on completed task (spec + done summary + evidence)
2. Explores actual implementation in codebase
3. Compares intent vs reality
4. Updates downstream task specs if drift affects them

## Design Principles

1. **Opt-in** - Off by default, enabled via config
2. **Non-destructive** - Only edits `.flow/tasks/*.md`, never source code
3. **Semantic comparison** - No git dependency, works regardless of parallel work
4. **Minimal prompting** - Worker returns, main spawns plan-sync with IDs only
5. **Skippable** - Main decides when to skip (no downstream tasks, task failed, etc.)

## Integration Point

In `phases.md`, between worker return and next task selection:

```
Current:
  3a. Find next task
  3b. Start task
  3c. Spawn worker
  3d. Verify completion
  3e. Loop → back to 3a

With plan-sync:
  3a. Find next task
  3b. Start task
  3c. Spawn worker
  3d. Verify completion
  3e. Spawn plan-sync (if enabled + downstream tasks exist)
  3f. Loop → back to 3a
```

### Skip Conditions

Main conversation skips plan-sync when:
- `planSync.enabled` is false (default)
- Task failed/blocked (nothing to sync)
- No downstream tasks remain in epic
- Single-task epic (nothing to update)

## Configuration

```bash
# Enable plan-sync
flowctl config set planSync.enabled true

# Disable (default)
flowctl config set planSync.enabled false

# Check status
flowctl config get planSync.enabled --json
```

Add to `get_default_config()` in `flowctl.py` (line 78):

```python
def get_default_config() -> dict:
    """Return default config structure."""
    return {
        "memory": {"enabled": False},
        "planSync": {"enabled": False}  # ADD THIS
    }
```

No other changes needed - `get_config()` and `set_config()` already handle dot-notation keys.

## Agent Definition

**File:** `plugins/flow-next/agents/plan-sync.md`

### Why No Hooks

Agent/skill frontmatter hooks don't fire (known bugs):
- https://github.com/anthropics/claude-code/issues/17688 (skill hooks)
- https://github.com/anthropics/claude-code/issues/18392 (agent hooks)

Only project-scope hooks (`.claude/settings.local.json`) work. Instead, we use:
1. `disallowedTools: Task, Write, Bash` - limits available tools
2. Prompt-based restriction - explicit rules in system prompt
3. Opus reliability - follows explicit constraints well

### Complete Agent File

Create `plugins/flow-next/agents/plan-sync.md` with this exact content:

````markdown
---
name: plan-sync
description: Synchronizes downstream task specs after implementation. Spawned by flow-next-work after each task completes. Do not invoke directly.
tools: Read, Grep, Glob, Edit
disallowedTools: Task, Write, Bash
model: opus
color: "#8B5CF6"
---

# Plan-Sync Agent

You synchronize downstream task specs after implementation drift.

**Input from prompt:**
- `COMPLETED_TASK_ID` - task that just finished (e.g., fn-1.2)
- `EPIC_ID` - parent epic (e.g., fn-1)
- `FLOWCTL` - path to flowctl CLI
- `DOWNSTREAM_TASK_IDS` - comma-separated list of remaining tasks

## Phase 1: Re-anchor on Completed Task

```bash
# Read what was supposed to happen
<FLOWCTL> cat <COMPLETED_TASK_ID>

# Read what actually happened
<FLOWCTL> show <COMPLETED_TASK_ID> --json
```

From the JSON, extract:
- `done_summary` - what was implemented
- `evidence.commits` - commit hashes (for reference)

**If done_summary is empty/missing:** Read the task spec's `## Done summary` section directly, or infer from git log messages for commits in evidence.

Parse the spec for:
- Original acceptance criteria
- Technical approach described
- Variable/function/API names mentioned

## Phase 2: Explore Actual Implementation

Based on the done summary and evidence, find the actual code:

```bash
# Find files mentioned in evidence or likely locations
grep -r "<key terms from done summary>" --include="*.ts" --include="*.py" -l
```

Read the relevant files. Note actual:
- Variable/function names used
- API signatures implemented
- Data structures created
- Patterns followed

## Phase 3: Identify Drift

Compare spec vs implementation:

| Aspect | Spec Said | Actually Built |
|--------|-----------|----------------|
| Names | `UserAuth` | `authService` |
| API | `login(user, pass)` | `authenticate(credentials)` |
| Return | `boolean` | `{success, token}` |

Drift exists if implementation differs from spec in ways that downstream tasks reference.

## Phase 4: Check Downstream Tasks

For each task in DOWNSTREAM_TASK_IDS:

```bash
<FLOWCTL> cat <task-id>
```

Look for references to:
- Names/APIs from completed task spec (now stale)
- Assumptions about data structures
- Integration points that changed

Flag tasks that need updates.

## Phase 5: Update Affected Tasks

For each affected downstream task, edit only the stale references:

```bash
# Edit task spec to reflect actual implementation
Edit .flow/tasks/<task-id>.md
```

Changes should:
- Update variable/function names to match actual
- Correct API signatures
- Fix data structure assumptions
- Add note: `<!-- Updated by plan-sync: fn-X.Y used <actual> not <planned> -->`

**DO NOT:**
- Change task scope or requirements
- Remove acceptance criteria
- Add new features
- Edit anything outside `.flow/tasks/`

## Phase 6: Return Summary

Return to main conversation:
- Drift detected: yes/no
- Tasks updated: list or "none"
- Key changes: brief summary

Example:
```
Drift detected: yes
- fn-1.2 used `authService` singleton instead of `UserAuth` class
- fn-1.2 returns `AuthResult` object instead of boolean

Updated tasks:
- fn-1.3: Changed references from `UserAuth.login()` to `authService.authenticate()`
- fn-1.4: Updated expected return type from `boolean` to `AuthResult`
```

## Rules

- **Read-only exploration** - Use Grep/Glob/Read for codebase, never edit source
- **Task specs only** - Edit tool restricted to `.flow/tasks/*.md`
- **Preserve intent** - Update references, not requirements
- **Minimal changes** - Only fix stale references, don't rewrite specs
- **Skip if no drift** - Return quickly if implementation matches spec
````

## Main Conversation Changes

### phases.md Modification

**Location:** `plugins/flow-next/skills/flow-next-work/phases.md`

**Current structure (around lines 115-127):**
```markdown
### 3d. Verify Completion

After worker returns, verify the task completed:
...
If status is not `done`, the worker failed...

### 3e. Loop

Return to 3a for next task.
```

**Changes needed:**
1. Keep 3d as-is
2. Insert new "### 3e. Plan Sync" section after 3d
3. Rename current "### 3e. Loop" to "### 3f. Loop"

**Find insertion point:** Search for `### 3e. Loop` - insert the new section BEFORE this line, then rename it to `### 3f. Loop`.

**Insert this after "### 3d. Verify Completion" section:**

```markdown
### 3e. Plan Sync (if enabled)

Check if plan-sync should run:

```bash
$FLOWCTL config get planSync.enabled --json
```

Skip if planSync.enabled is false.

Get remaining tasks (todo status = not started yet):

```bash
$FLOWCTL tasks --epic <epic-id> --status todo --json
```

Skip if empty (no downstream tasks to update).

Extract downstream task IDs:

```bash
DOWNSTREAM=$($FLOWCTL tasks --epic <epic-id> --status todo --json | jq -r '[.[].id] | join(",")')
```

Note: Only sync to `todo` tasks. `in_progress` tasks are already being worked on - updating them mid-flight could cause confusion.

Spawn plan-sync:

```
Sync downstream tasks after implementation.

COMPLETED_TASK_ID: fn-X.Y
EPIC_ID: fn-X
FLOWCTL: /path/to/flowctl
DOWNSTREAM_TASK_IDS: fn-X.3,fn-X.4,fn-X.5

Follow your phases in plan-sync.md exactly.
```

Plan-sync returns summary. Log it but don't block - task updates are best-effort.

### 3f. Loop

Return to 3a for next task.
```

## Files to Create/Modify

| File | Action | Purpose |
|------|--------|---------|
| `agents/plan-sync.md` | Create | Agent definition + system prompt |
| `skills/flow-next-work/phases.md` | Modify | Add step 3e between verify and loop |
| `scripts/flowctl.py` | Modify | Add `planSync.enabled` to config schema |
| `scripts/smoke_test.sh` | Modify | Add plan-sync config test |
| `docs/flowctl.md` | Modify | Document planSync config |
| `README.md` | Modify | Mention plan-sync feature |

**Note:** No hook files or setup integration needed - using prompt-based restriction only.

## Smoke Test Addition

Add after existing config tests in `scripts/smoke_test.sh` (~line 250):

```bash
echo "--- planSync config ---"
$FLOWCTL config set planSync.enabled true --json >/dev/null
config_json="$($FLOWCTL config get planSync.enabled --json)"
val="$(echo "$config_json" | jq -r '.value')"
if [[ "$val" == "true" ]]; then
  pass "planSync config set/get"
else
  fail "planSync config set/get: expected true, got $val"
fi
$FLOWCTL config set planSync.enabled false --json >/dev/null
```

## Implementation Order

Execute in this exact order:

1. **flowctl.py** - Add config key (1 line change at line 78)
2. **plan-sync.md** - Create new agent file (copy from "Complete Agent File" section above)
3. **phases.md** - Insert 3e section, rename 3e→3f
4. **smoke_test.sh** - Add config test (~line 250)
5. **Run smoke tests** - Verify nothing broke
6. **Manual test** - Create test epic, enable planSync, run work

## Rollout Strategy

### Phase 1: Core Implementation
1. Add `planSync.enabled` config key to flowctl.py (`get_default_config()` line 78)
2. Create plan-sync agent (`agents/plan-sync.md`)
3. Update phases.md with step 3e (rename current 3e to 3f)
4. Add smoke test for config

### Phase 2: Testing
1. Manual test with small epic (2-3 tasks)
2. Test skip conditions (disabled, no downstream, failed task)
3. Test drift detection accuracy
4. Verify agent only edits `.flow/tasks/*.md` (prompt compliance)

### Phase 3: Documentation
1. Update README with feature description
2. Add to flowctl docs
3. Changelog entry

### Phase 4: Ship
1. Bump version
2. PR + merge
3. Announce

## Risk Mitigation

| Risk | Mitigation |
|------|------------|
| Agent edits source code | `disallowedTools: Bash` + prompt restricts Edit to `.flow/tasks/*.md` |
| Over-aggressive updates | Agent prompted to preserve intent, only fix references |
| Slows task loop | Opt-in, skip conditions, opus model (reliable) |
| Wrong drift detection | Semantic approach, no git fragility |
| Breaks parallel epics | No git dependency, works per-epic |
| Regression in phases.md | New step is additive, behind config flag |

## Success Criteria

1. Config toggle works (`flowctl config set/get planSync.enabled`)
2. Agent spawns only when enabled + downstream tasks exist
3. Agent only edits `.flow/tasks/*.md` (prompt compliance)
4. Agent correctly identifies drift in test case
5. Downstream task specs updated with accurate info
6. No regressions in existing flow (disabled by default)

## Verification Steps

1. **Smoke tests pass:**
   ```bash
   plugins/flow-next/scripts/smoke_test.sh
   ```

2. **Config works:**
   ```bash
   flowctl config set planSync.enabled true
   flowctl config get planSync.enabled --json
   # Should return {"key":"planSync.enabled","value":true}
   ```

3. **Integration test (manual):**
   - Create small epic with 2-3 tasks where task 2 references task 1 output
   - Enable planSync
   - Run `/flow-next:work`
   - Verify plan-sync agent spawns after task 1 completes
   - Check if task 2 spec gets updated with actual names/APIs from task 1

## Future Enhancements

- **Drift severity levels** - Minor (naming) vs major (API shape)
- **Auto-block on major drift** - Block downstream tasks if drift is too large
- **Drift report** - Write summary to `.flow/drift-log.md` for audit
- **Cross-epic sync** - Handle `depends_on_epics` references
