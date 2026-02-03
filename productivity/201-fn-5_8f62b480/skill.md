# fn-5 TUI: ETA calculation

## Overview
Estimate time remaining for Ralph run based on task completion rate. Display in header or output panel.

## Scope
- Calculate ETA from tasks completed vs total
- Display formatted duration (e.g., "~8m", "~1h 20m")
- Update as tasks complete

## Approach
```typescript
function calculateETA(elapsed: number, completed: number, total: number): string {
  if (completed === 0) return '...';
  const avgPerTask = elapsed / completed;
  const remaining = total - completed;
  const etaSeconds = avgPerTask * remaining;
  return formatDuration(etaSeconds);
}
```

Display options:
- Header: `fn-1.3 「Add validation」 ETA ~8m  14:32`
- Output panel header: `─Iteration 3──────────────────── ETA ~8m ─`

## Quick commands
- `bun test`

## Acceptance
- [ ] ETA displays in UI when tasks > 0 completed
- [ ] ETA updates when task completes
- [ ] Shows "..." when no tasks completed yet
- [ ] Reasonable accuracy (within 20% of actual)

## References
- flow-next-tui spec: `plans/flow-next-tui-spec.md`
- Depends on MVP timer and task status tracking
