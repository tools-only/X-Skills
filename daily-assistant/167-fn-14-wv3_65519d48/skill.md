# TUI: Ralph Control Integration

## Overview

Integrate Ralph async control commands into flow-next-tui, allowing users to monitor and control Ralph runs from the terminal UI.

**Depends on:** fn-13-pxj (Ralph Async Control)

## Scope

**In scope:**
- Display pause/stop state in run list and detail views
- Show `completion_reason` for finished runs
- Add keybindings for pause/resume/stop actions
- Real-time state updates (poll or watch sentinel files)

**Out of scope:**
- Task reset from TUI (use CLI for now)
- Multiple concurrent runs

## Key Integration Points

From research (flow-next-tui codebase):
- `src/lib/runs.ts` — run discovery, `isRunActive()` at L162-185
- `src/lib/types.ts` — `TaskStatus`, need new `RunState` values
- Need to detect PAUSE/STOP sentinel files in run directories

## Approach

1. **Extend run state detection:**
   - Add `getRunControlState(runPath)` → `{paused: boolean, stopped: boolean, completion_reason?: string}`
   - Update `isRunActive()` to account for completion marker

2. **UI changes:**
   - Run list: show [PAUSED] / [STOPPED] badges
   - Run detail: show control state + completion reason
   - Add keybindings: `p` pause, `r` resume, `s` stop

3. **Actions:**
   - Call `flowctl ralph pause/resume/stop` via subprocess
   - Refresh state after action

## Quick Commands

```bash
# Dev server
cd flow-next-tui && bun dev

# Build
cd flow-next-tui && bun run build
```

## Acceptance

- [ ] Run list shows pause/stop state badges
- [ ] Run detail shows completion_reason for finished runs
- [ ] Keybindings work for pause/resume/stop
- [ ] State updates after control actions
- [ ] Works with sentinel files from fn-13-pxj

## References

- `flow-next-tui/src/lib/runs.ts` — run discovery
- `flow-next-tui/src/lib/types.ts` — type definitions
- fn-13-pxj spec — sentinel file locations and semantics
