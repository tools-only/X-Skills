# fn-13-pxj.2 Add sentinel file checks to ralph.sh

## Description
Add sentinel file checking to ralph.sh template so Ralph can be paused/stopped externally.

**CRITICAL**: Also ensure ALL exit paths write completion marker to progress.txt. Currently, some paths (e.g., `status == "none"`) exit without writing completion marker, making finished runs look "active" forever.

### Completion Marker Format

**CANONICAL MARKER**: `promise=COMPLETE` (plain text, NOT XML tagged)

This must match what `find_active_runs()` searches for in task fn-13-pxj.1. The detection checks for substring `promise=COMPLETE`, so the marker must be exactly that.

### Implementation

**In `plugins/flow-next/skills/flow-next-ralph-init/templates/ralph.sh`:**

1. Add `write_completion_marker()` function:
   ```bash
   write_completion_marker() {
     local reason="${1:-DONE}"
     echo "" >> "$PROGRESS_FILE"
     echo "completion_reason=$reason" >> "$PROGRESS_FILE"
     echo "promise=COMPLETE" >> "$PROGRESS_FILE"  # MUST match find_active_runs() detection
   }
   ```

2. Add `check_sentinels()` function:
   ```bash
   check_sentinels() {
     local pause_file="$RUN_DIR/PAUSE"
     local stop_file="$RUN_DIR/STOP"

     # Check for stop (do NOT remove file - keep for audit)
     if [[ -f "$stop_file" ]]; then
       log "STOP sentinel detected, exiting gracefully"
       write_completion_marker "STOPPED"
       exit 0
     fi

     # Check for pause (log once, not every loop)
     if [[ -f "$pause_file" ]]; then
       log "PAUSED - waiting for resume..."
       while [[ -f "$pause_file" ]]; do
         sleep 5
       done
       log "Resumed"
     fi
   }
   ```

3. Update ALL exit paths to call `write_completion_marker()`:
   - `status == "none"` (no work): `write_completion_marker "NO_WORK"`
   - Max iterations reached: `write_completion_marker "MAX_ITERATIONS"` (before exit, even if non-zero)
   - Normal completion: `write_completion_marker "DONE"`
   - FAIL paths: `write_completion_marker "FAILED"` (before exit 1)

4. Call `check_sentinels()` in main loop:
   - At start of each iteration (before work selection)
   - After Claude returns (before sleep/next iteration)

### Sentinel Semantics
- File presence = state active
- PAUSE: log once on enter, loop with 5s sleep until removed, log "Resumed" on exit
- STOP: log message, write completion marker, exit 0 (do NOT remove STOP file)
- Checking happens at iteration boundary only (won't interrupt running Claude)

### Key Files
- `plugins/flow-next/skills/flow-next-ralph-init/templates/ralph.sh`
## Acceptance
- [ ] `write_completion_marker()` function added
- [ ] `check_sentinels()` function added to ralph.sh template
- [ ] PAUSE sentinel causes loop to sleep, logs once on enter + once on resume
- [ ] STOP sentinel writes completion marker then exits (file NOT removed)
- [ ] Function called at start of each iteration (before work selection)
- [ ] Function called after Claude returns (before sleep/next iteration)
- [ ] ALL exit paths call `write_completion_marker()` with appropriate reason
- [ ] `status == "none"` exit writes "NO_WORK" marker
- [ ] Max iterations exit writes "MAX_ITERATIONS" marker
## Done summary
- Added `write_completion_marker()` - writes `promise=COMPLETE` + reason to progress.txt
- Added `check_sentinels()` - PAUSE loops with 5s sleep, STOP exits gracefully
- Updated 4 exit paths: NO_WORK, DONE, FAILED, MAX_ITERATIONS
- Sentinel checks at iteration start + after Claude returns

Why: Enable external pause/stop control + proper active run detection

Verification: CI (31/31), smoke (28/28), bash syntax OK
## Evidence
- Commits: a7a9720c69c27c323f52eadce453ad38e275289a
- Tests: plugins/flow-next/scripts/ci_test.sh, plugins/flow-next/scripts/smoke_test.sh
- PRs: