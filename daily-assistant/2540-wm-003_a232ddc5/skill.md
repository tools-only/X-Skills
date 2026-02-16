# WM-003: Extend StateManager with metrics methods

## Overview
Add new methods to StateManager for recording timestamps and storing/retrieving metrics.

## Files to Modify
- `zerg/state.py`

## Requirements

### New Imports
At the top of the file, add to the imports from zerg.types:
```python
from zerg.types import ExecutionEvent, WorkerState, FeatureMetrics
```

### New Methods

Add these methods to the StateManager class:

```python
def record_task_claimed(self, task_id: str, worker_id: int) -> None:
    """Record when a task was claimed by a worker.

    Sets the claimed_at timestamp for task metrics tracking.

    Args:
        task_id: Task identifier
        worker_id: Worker that claimed the task
    """
    with self._lock:
        if "tasks" not in self._state:
            self._state["tasks"] = {}

        if task_id not in self._state["tasks"]:
            self._state["tasks"][task_id] = {}

        self._state["tasks"][task_id]["claimed_at"] = datetime.now().isoformat()
        self._state["tasks"][task_id]["worker_id"] = worker_id

    self.save()
    logger.debug(f"Task {task_id} claimed by worker {worker_id}")


def record_task_duration(self, task_id: str, duration_ms: int) -> None:
    """Record task execution duration.

    Args:
        task_id: Task identifier
        duration_ms: Execution duration in milliseconds
    """
    with self._lock:
        if task_id in self._state.get("tasks", {}):
            self._state["tasks"][task_id]["duration_ms"] = duration_ms

    self.save()
    logger.debug(f"Task {task_id} duration: {duration_ms}ms")


def store_metrics(self, metrics: "FeatureMetrics") -> None:
    """Store computed metrics to state.

    Args:
        metrics: FeatureMetrics to persist
    """
    with self._lock:
        self._state["metrics"] = metrics.to_dict()

    self.save()
    logger.debug("Stored feature metrics")


def get_metrics(self) -> "FeatureMetrics | None":
    """Retrieve stored metrics.

    Returns:
        FeatureMetrics if available, None otherwise
    """
    with self._lock:
        metrics_data = self._state.get("metrics")
        if not metrics_data:
            return None

        # Import here to avoid circular import
        from zerg.types import FeatureMetrics
        return FeatureMetrics.from_dict(metrics_data)
```

### Update set_worker_ready

The existing `set_worker_ready` method already sets `ready_at`. Verify it looks like this:

```python
def set_worker_ready(self, worker_id: int) -> None:
    """Mark a worker as ready to receive tasks.

    Args:
        worker_id: Worker identifier
    """
    with self._lock:
        worker_data = self._state.get("workers", {}).get(str(worker_id), {})
        if worker_data:
            worker_data["status"] = WorkerStatus.READY.value
            worker_data["ready_at"] = datetime.now().isoformat()

    self.save()
    logger.debug(f"Worker {worker_id} marked ready")
```

### Update claim_task

Modify the existing `claim_task` method to also record claimed_at:

```python
def claim_task(self, task_id: str, worker_id: int) -> bool:
    """Attempt to claim a task for a worker.

    Uses file-based locking to prevent race conditions.

    Args:
        task_id: Task to claim
        worker_id: Worker claiming the task

    Returns:
        True if claim succeeded
    """
    with self._lock:
        # Reload state to get latest
        self.load()

        task_state = self._state.get("tasks", {}).get(task_id, {})
        current_status = task_state.get("status", TaskStatus.PENDING.value)

        # Can only claim pending tasks
        if current_status not in (TaskStatus.TODO.value, TaskStatus.PENDING.value):
            return False

        # Check if already claimed by a different worker
        existing_worker = task_state.get("worker_id")
        if existing_worker is not None and existing_worker != worker_id:
            return False

        # Claim it - now also records claimed_at
        self.set_task_status(task_id, TaskStatus.CLAIMED, worker_id=worker_id)
        self.record_task_claimed(task_id, worker_id)  # Add this line
        logger.info(f"Worker {worker_id} claimed task {task_id}")
        return True
```

### Update _create_initial_state

Add metrics to initial state:

```python
def _create_initial_state(self) -> dict[str, Any]:
    """Create initial state structure.

    Returns:
        Initial state dictionary
    """
    return {
        "feature": self.feature,
        "started_at": datetime.now().isoformat(),
        "current_level": 0,
        "tasks": {},
        "workers": {},
        "levels": {},
        "execution_log": [],
        "metrics": None,  # Add this line
        "paused": False,
        "error": None,
    }
```

## Verification

```bash
python -c "from zerg.state import StateManager; s = StateManager('test-metrics-check'); print('StateManager OK')"
```

## Acceptance Criteria
- [ ] record_task_claimed method sets claimed_at timestamp
- [ ] record_task_duration method sets duration_ms field
- [ ] store_metrics method persists FeatureMetrics to state
- [ ] get_metrics method returns FeatureMetrics or None
- [ ] claim_task now records claimed_at timestamp
- [ ] Initial state includes metrics field
- [ ] All methods are thread-safe (use self._lock)
