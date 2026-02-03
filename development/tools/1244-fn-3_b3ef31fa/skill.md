# fn-3 Async Mode: CLI Control for External Agents

## Overview

Ralph runs as a foreground process. Users want to:
- Start Ralph, walk away, monitor from phone/Discord via Clawdbot
- Pause/resume without killing the process
- Rollback failed tasks and retry without restarting everything
- Check status from any CLI tool

## Scope

Add CLI commands and sentinel file mechanism. Any tool that can run shell commands (Clawdbot, cron, future TUI) can control Ralph.

**Not in scope:** Session attachment (use tmux), web UI, TUI. Keep it simple and composable.

## Approach

### 1. Task Reset Command

```bash
flowctl task reset fn-1.3              # done → pending
flowctl task reset --cascade fn-1.3    # + reset all dependents
```

- Validates task exists
- Only resets if status is `done` or `blocked`
- `--cascade` finds all tasks with `depends_on` containing this task (recursively)
- Updates `updated_at` timestamp

### 2. Status Command

```bash
flowctl status                         # Human-readable summary
flowctl status --json                  # Machine-readable for bots
```

Output includes:
- Active epic(s) with progress (X/Y tasks done)
- Current task (if any in_progress)
- Last completed task + timestamp
- Any blocked tasks with reason
- Ralph state (running/paused/stopped) if sentinel files present

### 3. Sentinel File Mechanism

Ralph checks between tasks for:
- `scripts/ralph/PAUSE` → finish current task, then wait
- `scripts/ralph/STOP` → finish current task, then exit cleanly

```bash
# Pause Ralph
touch scripts/ralph/PAUSE

# Resume
rm scripts/ralph/PAUSE

# Stop gracefully
touch scripts/ralph/STOP
```

Ralph loop changes:
```bash
# After each task, before next:
if [[ -f scripts/ralph/STOP ]]; then
  log "STOP requested, exiting"
  rm scripts/ralph/STOP
  exit 0
fi

while [[ -f scripts/ralph/PAUSE ]]; then
  log "PAUSED - waiting for PAUSE removal"
  sleep 5
done
```

### 4. Convenience wrappers

```bash
flowctl ralph pause    # touch PAUSE
flowctl ralph resume   # rm PAUSE
flowctl ralph stop     # touch STOP
```

## Quick commands

- `plugins/flow-next/scripts/smoke_test.sh`
- `plugins/flow-next/scripts/ralph_smoke_test.sh`

## Acceptance

- [ ] `flowctl task reset fn-X.Y` marks task pending
- [ ] `flowctl task reset --cascade fn-X.Y` resets task + dependents
- [ ] `flowctl status` shows current state
- [ ] `flowctl status --json` outputs machine-readable format
- [ ] `touch scripts/ralph/PAUSE` pauses Ralph after current task
- [ ] `rm scripts/ralph/PAUSE` resumes Ralph
- [ ] `touch scripts/ralph/STOP` stops Ralph gracefully
- [ ] Smoke tests cover new commands

## References

- Twitter thread with @tiagoefreitas on async monitoring
- Clawdbot integration use case
