# fn-13-pxj.3 Add flowctl ralph pause/resume/stop/status

## Description
Add `flowctl ralph` subcommand group with pause/resume/stop/status.

### Implementation

**In flowctl.py:**

1. Add `find_active_run()` helper (similar to find_active_runs but returns single run):
   ```python
   def find_active_run(run_id: str | None = None, use_json: bool = False) -> Path:
       """Find active run directory. Auto-detect if run_id is None."""
       runs = find_active_runs()
       if run_id:
           matches = [r for r in runs if r[0] == run_id]
           if not matches:
               error_exit(f"Run {run_id} not found or not active", use_json)
           return matches[0][1]  # run_dir
       if len(runs) == 0:
           error_exit("No active runs", use_json)
       if len(runs) > 1:
           ids = ', '.join(r[0] for r in runs)
           error_exit(f"Multiple active runs, specify --run: {ids}", use_json)
       return runs[0][1]
   ```

2. Add command handlers (using `error_exit()` not `fail()`):
   ```python
   def cmd_ralph_pause(args):
       run_dir = find_active_run(args.run, args.json)
       (run_dir / "PAUSE").touch()
       json_output({"success": True, "run": run_dir.name, "action": "paused"}) if args.json else print(f"Paused {run_dir.name}")

   def cmd_ralph_resume(args):
       run_dir = find_active_run(args.run, args.json)
       (run_dir / "PAUSE").unlink(missing_ok=True)
       json_output({"success": True, "run": run_dir.name, "action": "resumed"}) if args.json else print(f"Resumed {run_dir.name}")

   def cmd_ralph_stop(args):
       run_dir = find_active_run(args.run, args.json)
       (run_dir / "STOP").touch()
       # Note: STOP file is NOT removed by ralph.sh - kept for audit
       json_output({"success": True, "run": run_dir.name, "action": "stop_requested"}) if args.json else print(f"Stop requested for {run_dir.name}")

   def cmd_ralph_status(args):
       run_dir = find_active_run(args.run, args.json)
       paused = (run_dir / "PAUSE").exists()
       stopped = (run_dir / "STOP").exists()
       # Read progress.txt for iteration/current task
       # Output JSON or human-readable
   ```

3. Register nested subparsers:
   ```python
   p_ralph = subparsers.add_parser("ralph", help="Ralph control commands")
   ralph_sub = p_ralph.add_subparsers(dest="ralph_cmd", required=True)
   
   p_pause = ralph_sub.add_parser("pause", help="Pause Ralph run")
   p_pause.add_argument("--run", help="Run ID (auto-detect if single)")
   p_pause.add_argument("--json", action="store_true")
   p_pause.set_defaults(func=cmd_ralph_pause)
   # ... similar for resume, stop, status
   ```

### Key Behaviors
- STOP file is NOT removed by ralph.sh (kept for audit/inspection)
- `ralph status` shows paused=true if PAUSE exists, stopped=true if STOP exists
- All error paths use `error_exit(message, use_json=args.json)`

### Key Files
- `plugins/flow-next/scripts/flowctl.py`
## Acceptance
- [ ] `flowctl ralph pause` creates PAUSE sentinel in run directory
- [ ] `flowctl ralph resume` removes PAUSE sentinel
- [ ] `flowctl ralph stop` creates STOP sentinel
- [ ] `flowctl ralph status` shows run state (active/paused/stopped)
- [ ] Auto-detects run_id when only one active run
- [ ] Errors with run list when multiple active runs
- [ ] `--run <id>` flag works to specify run
- [ ] `--json` flag works for all subcommands
- [ ] Pause already-paused is no-op (success)
- [ ] Resume non-paused is no-op (success)
## Done summary
- Added `find_active_run()` helper for single-run auto-detection
- Added `cmd_ralph_pause/resume/stop/status` handlers
- Registered `flowctl ralph` subparser group

Why: Enable external control of Ralph runs via CLI

Verification: Manual tests for pause/resume/stop/status, auto-detect, error cases; CI 31/31
## Evidence
- Commits: 784fb1e828e3accc895356e4ec14476d465a42f1
- Tests: plugins/flow-next/scripts/ci_test.sh
- PRs: