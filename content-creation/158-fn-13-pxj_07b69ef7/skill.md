# Ralph Async Control

**Issue:** https://github.com/gmickel/gmickel-claude-marketplace/issues/14

## Overview

Add CLI commands to pause/resume/stop Ralph runs and reset tasks, enabling external agents (Clawdbot, GitHub Actions, etc.) to control runs without killing processes.

## Scope

**In scope:**
- `flowctl status [--json]` - show epic/task counts + active run info
- Sentinel files (PAUSE/STOP) in run directory
- `flowctl ralph pause/resume/stop/status [--run <run_id>]`
- `flowctl task reset <task_id> [--cascade]`
- `flowctl epic add-dep/rm-dep` - manage epic-level dependencies via CLI
- ralph.sh sentinel checking at iteration boundary
- ralph.sh completion marker on ALL exit paths (fix for active run detection)
- CI tests for all commands (jq-free, Python-based)
- Documentation updates (README, docs/flowctl.md)

**Out of scope:**
- Session attachment (use tmux/screen)
- Web UI / TUI controls (TUI is separate package)
- Notifications (external agents poll status)
- Multiple concurrent runs

## Approach

### Sentinel File Mechanism
- Location: `scripts/ralph/runs/<run_id>/PAUSE` and `STOP`
- Semantics: file presence = state active; remove file to clear state
- Timing: ralph.sh checks at iteration boundary (after Claude call completes)
- No file locking needed (atomic touch/rm operations)

### Active Run Detection (CRITICAL)
Current ralph.sh exits without writing completion marker on some paths (e.g., `status == "none"`). Must fix ALL exit paths to write completion marker.

**Canonical marker string**: `promise=COMPLETE` (plain text, NOT XML)
- Written to progress.txt by `write_completion_marker()` in ralph.sh
- Detected by `find_active_runs()` in flowctl.py via substring search
- MUST be consistent between writer and reader

Exit reasons appended alongside marker:
- `completion_reason=DONE` - normal completion
- `completion_reason=NO_WORK` - no tasks available
- `completion_reason=STOPPED` - STOP sentinel detected
- `completion_reason=MAX_ITERATIONS` - iteration limit reached
- `completion_reason=FAILED` - error exit

### flowctl Commands

**`flowctl status [--json]`**
- Reads .flow/ state (epic/task counts)
- Scans `scripts/ralph/runs/*/progress.txt` for active runs
- Active = progress.txt exists AND lacks substring `promise=COMPLETE`
- Shows: active run_id, iteration, current epic/task, paused/stopped state

**`flowctl ralph <cmd> [--run <run_id>]`**
- `pause`: touch PAUSE sentinel
- `resume`: rm PAUSE sentinel
- `stop`: touch STOP sentinel (NOT removed on ralph exit; kept for audit)
- `status`: show run state (active/paused/stopped, iteration, current work)
- Auto-detect run_id if only one active; error with list if multiple
- All commands support `--json` flag with `error_exit(msg, use_json=args.json)` pattern

**`flowctl task reset <task_id> [--cascade]`**
- Sets status → `todo`
- Clears: `blocked_reason` (if present), `completed_at` (if present)
- Clears claim fields: `assignee`, `claimed_at`, `claim_note`
- Clears `evidence` from JSON + clears `## Evidence` section contents in spec (keeps heading)
- Cannot reset `in_progress` tasks (error)
- Cannot reset tasks in closed epic (error)
- `--cascade`: also reset dependent tasks (same epic only, skip in_progress)

### ralph.sh Changes
- Add `write_completion_marker()` function (writes `promise=COMPLETE`)
- Add `check_sentinels()` function
- Call after each iteration (after Claude exits, before sleep/next)
- PAUSE: loop with 5s sleep until PAUSE removed, log once on enter + on resume (not every loop)
- STOP: write completion marker to progress.txt, exit 0 (do NOT remove STOP file)
- **CRITICAL**: Ensure ALL exit paths write `promise=COMPLETE` to progress.txt

## Quick Commands

```bash
# Full CI test suite (includes new async control tests)
cd plugins/flow-next && ./scripts/ci_test.sh

# Manual smoke test: create mock run and test controls
mkdir -p scripts/ralph/runs/test-run
echo "iteration: 1" > scripts/ralph/runs/test-run/progress.txt
flowctl ralph pause --run test-run
ls scripts/ralph/runs/test-run/  # Should show PAUSE file
flowctl ralph resume --run test-run
rm -rf scripts/ralph/runs/test-run
```

## Acceptance

- [ ] `flowctl status` shows epic/task counts
- [ ] `flowctl status --json` returns valid JSON with run info
- [ ] `flowctl ralph pause/resume/stop/status` work with auto-detection
- [ ] ralph.sh respects PAUSE sentinel (pauses loop, logs once)
- [ ] ralph.sh respects STOP sentinel (writes completion marker, exits)
- [ ] ralph.sh writes `promise=COMPLETE` (plain text) on ALL exit paths
- [ ] `find_active_runs()` correctly detects active vs completed runs
- [ ] `flowctl task reset` changes done/blocked → todo
- [ ] `flowctl task reset` clears claim fields (assignee, claimed_at, claim_note)
- [ ] `flowctl task reset` clears evidence from JSON + spec (keeps heading)
- [ ] `flowctl task reset --cascade` resets dependent tasks
- [ ] `flowctl task reset` errors on in_progress tasks
- [ ] `flowctl epic add-dep/rm-dep` manage epic dependencies
- [ ] All CI tests pass (jq-free)
- [ ] Documentation updated (README, docs/flowctl.md)

## References

### Code Locations
- `plugins/flow-next/scripts/flowctl.py` - CLI implementation (argparse at L3711+)
- `plugins/flow-next/skills/flow-next-ralph-init/templates/ralph.sh` - main loop L725-963, run dir L488-501
- `plugins/flow-next/scripts/ci_test.sh` - test patterns (Python JSON parsing, no jq)

### Patterns
- Subparsers: `subparsers.add_parser()` + `set_defaults(func=)`
- Run active check: `progress.txt` lacking substring `promise=COMPLETE`
- Test helpers: `pass()` / `fail()` with color output
- Error handling: `error_exit(message, use_json=args.json)`

### Run Directory Structure
```
scripts/ralph/runs/<run_id>/
├── progress.txt      # Ralph progress log (MUST contain promise=COMPLETE when done)
├── attempts.json     # Task attempt counts
├── branches.json     # Branch tracking
├── run.json          # Scoped epics (if set)
├── receipts/         # Review receipts
├── PAUSE             # NEW: pause sentinel (removed on resume)
└── STOP              # NEW: stop sentinel (kept after exit for audit)
```
