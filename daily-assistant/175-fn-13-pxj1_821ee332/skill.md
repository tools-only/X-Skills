# fn-13-pxj.1 Add flowctl status command

## Description
Add `flowctl status [--json]` command to show epic/task counts plus active Ralph run info.

### Prerequisites
Task fn-13-pxj.2 must ensure ralph.sh writes `promise=COMPLETE` on ALL exit paths, otherwise status will incorrectly show finished runs as "active".

### Implementation

**In flowctl.py:**

1. Add `find_active_runs()` helper:
   - Scan `scripts/ralph/runs/*/progress.txt`
   - Run is active if: progress.txt exists AND does NOT contain "promise=COMPLETE"
   - Return list of (run_id, run_dir, progress_info)

2. Add `cmd_status()`:
   - Get epic/task counts via existing functions
   - Call `find_active_runs()` for Ralph state
   - Check for PAUSE/STOP sentinel files in active runs
   - Output human-readable or JSON based on `--json`

3. Register parser:
   ```python
   p_status = subparsers.add_parser("status", help="Show .flow state and active runs")
   p_status.add_argument("--json", action="store_true")
   p_status.set_defaults(func=cmd_status)
   ```

### Output Format

**Human-readable:**
```
Epics: 2 open, 1 done
Tasks: 5 todo, 1 in_progress, 3 done, 1 blocked

Active runs:
  run-abc123 (iteration 5, working on fn-1-xyz.2)
```

**JSON:**
```json
{
  "epics": {"open": 2, "done": 1},
  "tasks": {"todo": 5, "in_progress": 1, "done": 3, "blocked": 1},
  "runs": [
    {
      "id": "run-abc123",
      "iteration": 5,
      "current_task": "fn-1-xyz.2",
      "paused": false,
      "stopped": false
    }
  ]
}
```

### Key Files
- `plugins/flow-next/scripts/flowctl.py` - add cmd_status and find_active_runs
## Acceptance
- [ ] `flowctl status` outputs epic/task counts
- [ ] `flowctl status` shows active Ralph runs (or "No active runs")
- [ ] `flowctl status --json` returns valid JSON
- [ ] JSON includes runs array with id, iteration, current_task, paused, stopped
- [ ] Works when no runs/ directory exists
- [ ] Works when runs/ is empty
## Done summary
- Added `find_active_runs()` helper to detect Ralph runs by scanning progress.txt
- Added `cmd_status()` showing epic/task counts + active runs with PAUSE/STOP state
- Registered `flowctl status [--json]` subparser

Why: Enable monitoring of .flow state and Ralph runs externally

Verification: CI tests pass (31/31), manual tests for active/paused/stopped states
## Evidence
- Commits: 4197e7e
- Tests: plugins/flow-next/scripts/ci_test.sh
- PRs: