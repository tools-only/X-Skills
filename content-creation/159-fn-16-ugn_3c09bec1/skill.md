# fn-16-ugn Manual Plan-Sync Trigger

## Overview

Add `/flow-next:sync` command to manually trigger plan-sync from a given task or epic, without requiring the full `/flow-next:work` loop.

Related: [Issue #43](https://github.com/gmickel/gmickel-claude-marketplace/issues/43)

## Scope

- New `/flow-next:sync` command (skill)
- Reuses existing `agents/plan-sync.md`
- Works with task ID or epic ID
- Supports `--dry-run` flag

## Design Notes

**Config behavior:**
- Manual command ignores `planSync.enabled` config - if user explicitly runs it, they want it
- Auto-trigger in `/flow-next:work` still respects the config

**Source task status:**
- Any status allowed (todo, in_progress, done, blocked)
- User may have edited spec/code manually - trust their intent
- Differs from auto-trigger which requires `done` status

**Safety:**
- Purely additive - doesn't modify existing skills/phases
- Agent works standalone (reads specs + git state, no work-loop context needed)
- Ralph unaffected - calls `/flow-next:work` which has its own flow

## Approach

**With task ID (fn-N.M):**
1. Validate task exists
2. Read the task's spec (evidence optional - agent compares spec vs code state)
3. Find all downstream `todo` or `blocked` tasks in same epic
4. If none: report "No downstream tasks to sync" and exit
5. Spawn plan-sync agent to update affected specs

**With epic ID (fn-N):**
1. Validate epic exists
2. Find all tasks in epic (any status)
3. Identify downstream `todo` or `blocked` tasks needing updates
4. If none: report "No downstream tasks to sync" and exit
5. Run plan-sync agent on the full set

**With --dry-run:**
- Agent shows proposed changes but doesn't write
- Output: "Would update: fn-N.M, fn-N.K" with diff preview

## Error Handling

Simple exit + guidance (no complex recovery):

| Case | Message |
|------|---------|
| No `.flow/` directory | "No .flow/ found. Run `flowctl init` first." |
| Task/epic not found | "Task fn-X.Y not found. Run `flowctl list` to see tasks." |
| No downstream tasks | "No downstream tasks to sync (all done or none exist)." |
| Invalid ID format | "Invalid ID format. Use fn-N (epic) or fn-N.M (task)." |

## Output

On success, show:
```
Plan-sync: fn-1.2 → downstream tasks

Scanned: 3 tasks (fn-1.3, fn-1.4, fn-1.5)
Drift detected: 2 tasks
Updated: fn-1.3, fn-1.5
No changes: fn-1.4

Summary: <agent's drift summary>
```

On dry-run:
```
Plan-sync: fn-1.2 → downstream tasks (DRY RUN)

Would update:
- fn-1.3: Update API endpoint reference (POST → PUT)
- fn-1.5: Update field name (userId → accountId)

No files modified.
```

## Quick commands

```bash
plugins/flow-next/scripts/smoke_test.sh
```

## Acceptance

- [ ] `/flow-next:sync fn-N.M` updates downstream tasks based on source task
- [ ] `/flow-next:sync fn-N` scans whole epic for drift
- [ ] `/flow-next:sync fn-N.M --dry-run` shows changes without writing
- [ ] Works with any source task status (not just done)
- [ ] Ignores `planSync.enabled` config (manual = always run)
- [ ] Includes `blocked` tasks in target set (not just todo)
- [ ] Clear error messages for missing .flow/, invalid ID, no tasks
- [ ] Reuses existing plan-sync agent (no duplication)
- [ ] Smoke test passes

## File Structure

```
plugins/flow-next/skills/flow-next-sync/
  SKILL.md       # Main skill (triggers on /flow-next:sync)
```

## References

- `agents/plan-sync.md` - existing agent to reuse
- `skills/flow-next-work/phases.md` section 3e - current auto-trigger logic
- Issue #43 - broader spec-sync feature request
