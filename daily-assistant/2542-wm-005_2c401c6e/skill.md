# WM-005: Integrate metrics into orchestrator

## Overview
Update the orchestrator to record timestamps at key lifecycle points and trigger metrics computation.

## Files to Modify
- `zerg/orchestrator.py`

## Requirements

### New Import
Add at the top with other imports:
```python
from zerg.metrics import MetricsCollector
```

### Record Worker Timestamps

In the `_spawn_worker` method (or wherever workers are spawned), ensure `started_at` is recorded:

```python
# When creating worker state
worker_state = WorkerState(
    worker_id=worker_id,
    status=WorkerStatus.INITIALIZING,
    started_at=datetime.now(),  # Ensure this is set
    # ... other fields
)
```

### Record Worker Ready

When a worker becomes ready, the orchestrator should call:
```python
self._state.set_worker_ready(worker_id)
```

This already records `ready_at` in the state.

### Record Task Claimed

When a task is claimed by a worker, ensure `record_task_claimed` is called. This may already happen in the state manager's `claim_task` method. Verify the orchestrator uses it:

```python
# In task assignment flow
if self._state.claim_task(task_id, worker_id):
    # Task claimed successfully
    # claimed_at is now recorded
    ...
```

### Compute Task Duration

When a task completes, compute and record duration:

```python
def _handle_task_completion(self, task_id: str, worker_id: int) -> None:
    """Handle successful task completion."""
    # Get task start time
    task_state = self._state._state.get("tasks", {}).get(task_id, {})
    started_at = task_state.get("started_at")

    if started_at:
        # Calculate duration
        from zerg.metrics import duration_ms
        now = datetime.now()
        task_duration = duration_ms(started_at, now)
        if task_duration:
            self._state.record_task_duration(task_id, task_duration)

    # Continue with existing completion logic
    self._state.set_task_status(task_id, TaskStatus.COMPLETE, worker_id=worker_id)
    # ... rest of completion handling
```

### Trigger Metrics Computation on Level Complete

In the level completion handler, compute and store metrics:

```python
def _on_level_complete_handler(self, level: int) -> None:
    """Handle level completion."""
    # Existing level completion logic...

    # Compute and store metrics
    collector = MetricsCollector(self._state)
    metrics = collector.compute_feature_metrics()
    self._state.store_metrics(metrics)

    # Log metrics summary
    logger.info(
        f"Level {level} complete: "
        f"{metrics.tasks_completed}/{metrics.tasks_total} tasks, "
        f"{metrics.total_duration_ms}ms total"
    )
```

### Add Metrics to Status Method

If the orchestrator has a `status()` or `get_status()` method, include metrics:

```python
def status(self) -> dict:
    """Get current execution status with metrics."""
    # Compute fresh metrics
    collector = MetricsCollector(self._state)
    metrics = collector.compute_feature_metrics()

    return {
        "feature": self.feature,
        "current_level": self._state.get_current_level(),
        "paused": self._state.is_paused(),
        "error": self._state.get_error(),
        "metrics": metrics.to_dict(),
        # ... other status fields
    }
```

## Verification

```bash
python -c "from zerg.orchestrator import Orchestrator; print('Orchestrator import OK')"
```

## Acceptance Criteria
- [ ] MetricsCollector is imported from zerg.metrics
- [ ] Worker started_at is recorded when spawned
- [ ] Worker ready_at is recorded when ready (via set_worker_ready)
- [ ] Task claimed_at is recorded when claimed
- [ ] Task duration_ms is computed and stored on completion
- [ ] MetricsCollector.compute_feature_metrics() called on level complete
- [ ] Metrics stored to state via store_metrics()
- [ ] Status output includes metrics
