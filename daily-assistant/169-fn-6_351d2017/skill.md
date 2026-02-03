# fn-6 TUI: Multi-run support

## Overview
Tab key cycles between Ralph runs in current project. Users may have multiple runs (interrupted, retried, different epics).

## Scope
- Scan all runs in `scripts/ralph/runs/`
- Tab cycles forward, Shift+Tab cycles backward
- Status bar shows run position (e.g., "2/5")
- State reloads on switch (tasks, output, scroll)

## Approach
1. On startup, scan `scripts/ralph/runs/` for all runs
2. Sort by date (newest first)
3. Store as `runs[]` array with `currentRunIndex`
4. Tab increments index (wraps around)
5. On switch: reload tasks, output buffer, reset scroll

## Quick commands
- `bun test`

## Acceptance
- [ ] Tab switches to next run
- [ ] Shift+Tab switches to previous run
- [ ] Status bar shows current position (2/5)
- [ ] Task list updates on switch
- [ ] Output buffer reloads on switch
- [ ] Wraps around at boundaries

## References
- flow-next-tui spec: `plans/flow-next-tui-spec.md`
- Depends on MVP run discovery
