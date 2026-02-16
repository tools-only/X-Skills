# DC-009: Wire ContainerLauncher to Orchestrator

**Level**: 4 | **Critical Path**: Yes ⭐ | **Estimate**: 35 min
**Dependencies**: DC-007, DC-008

## Objective

Complete the orchestrator integration so that `_spawn_worker`, `_poll_workers`, and `_terminate_worker` properly use `ContainerLauncher` when in container mode.

## Files Owned

- `zerg/orchestrator.py` (modify)

## Files to Read

- `zerg/launcher.py` (ContainerLauncher methods)

## Implementation

### 1. Update _spawn_worker()

The existing code already has some container logic. Update to use ContainerLauncher properly:

```python
def _spawn_worker(self, worker_id: int) -> WorkerState:
    """Spawn a single worker.

    Args:
        worker_id: Worker identifier

    Returns:
        WorkerState for the spawned worker
    """
    logger.info(f"Spawning worker {worker_id}")

    # Allocate port
    port = self.ports.allocate_one()

    # Create worktree
    wt_info = self.worktrees.create(self.feature, worker_id)

    # Check launcher type
    if isinstance(self.launcher, ContainerLauncher):
        # Container mode: spawn container and exec claude
        result = self.launcher.spawn_and_exec(
            worker_id=worker_id,
            feature=self.feature,
            worktree_path=wt_info.path,
            branch=wt_info.branch,
        )

        if not result.success:
            raise RuntimeError(f"Failed to spawn container worker: {result.error}")

        worker_state = WorkerState(
            worker_id=worker_id,
            status=WorkerStatus.RUNNING,
            port=port,
            container_id=result.handle.container_id if result.handle else None,
            worktree_path=str(wt_info.path),
            branch=wt_info.branch,
            started_at=datetime.now(),
        )
    else:
        # Subprocess mode: existing logic
        result = self.launcher.spawn(
            worker_id=worker_id,
            feature=self.feature,
            worktree_path=wt_info.path,
            branch=wt_info.branch,
        )

        if not result.success:
            raise RuntimeError(f"Failed to spawn worker: {result.error}")

        worker_state = WorkerState(
            worker_id=worker_id,
            status=WorkerStatus.RUNNING,
            port=port,
            worktree_path=str(wt_info.path),
            branch=wt_info.branch,
            started_at=datetime.now(),
        )

    self._workers[worker_id] = worker_state
    self.state.set_worker_state(worker_state)
    self.state.append_event("worker_started", {
        "worker_id": worker_id,
        "port": port,
        "mode": "container" if isinstance(self.launcher, ContainerLauncher) else "subprocess",
    })

    return worker_state
```

### 2. Update _poll_workers()

Already uses launcher.monitor() which works for both. Just ensure status mapping is correct:

```python
def _poll_workers(self) -> None:
    """Poll worker status and handle completions."""
    for worker_id, worker in list(self._workers.items()):
        # Monitor via launcher (works for both subprocess and container)
        status = self.launcher.monitor(worker_id)

        if status != worker.status:
            logger.debug(f"Worker {worker_id} status: {worker.status.value} -> {status.value}")

        if status == WorkerStatus.CRASHED:
            logger.error(f"Worker {worker_id} crashed")
            worker.status = WorkerStatus.CRASHED
            self.state.set_worker_state(worker)

            if worker.current_task:
                self._handle_task_failure(
                    worker.current_task,
                    worker_id,
                    "Worker crashed",
                )

        elif status == WorkerStatus.CHECKPOINTING:
            logger.info(f"Worker {worker_id} checkpointing")
            worker.status = WorkerStatus.CHECKPOINTING
            self.state.set_worker_state(worker)
            self._handle_worker_exit(worker_id)

        elif status == WorkerStatus.STOPPED:
            self._handle_worker_exit(worker_id)

        # Update health check
        worker.health_check_at = datetime.now()
```

### 3. Update _terminate_worker()

Use launcher.terminate() which works for both:

```python
def _terminate_worker(self, worker_id: int, force: bool = False) -> None:
    """Terminate a worker.

    Args:
        worker_id: Worker identifier
        force: Force termination
    """
    worker = self._workers.get(worker_id)
    if not worker:
        return

    logger.info(f"Terminating worker {worker_id} (mode: {self.get_launcher_mode()})")

    # Terminate via launcher (works for both modes)
    self.launcher.terminate(worker_id, force=force)

    # Delete worktree
    try:
        wt_path = self.worktrees.get_worktree_path(self.feature, worker_id)
        self.worktrees.delete(wt_path, force=True)
    except Exception as e:
        logger.warning(f"Failed to delete worktree for worker {worker_id}: {e}")

    # Release port
    if worker.port:
        self.ports.release(worker.port)

    # Update state
    worker.status = WorkerStatus.STOPPED
    self.state.set_worker_state(worker)
    self.state.append_event("worker_stopped", {
        "worker_id": worker_id,
        "mode": self.get_launcher_mode(),
    })

    del self._workers[worker_id]
```

### 4. Add container-specific status method

```python
def status(self) -> dict[str, Any]:
    """Get current orchestration status.

    Returns:
        Status dictionary
    """
    level_status = self.levels.get_status()

    return {
        "feature": self.feature,
        "running": self._running,
        "launcher_mode": self.get_launcher_mode(),  # NEW
        "current_level": level_status["current_level"],
        "progress": {
            "total": level_status["total_tasks"],
            "completed": level_status["completed_tasks"],
            "failed": level_status["failed_tasks"],
            "in_progress": level_status["in_progress_tasks"],
            "percent": level_status["progress_percent"],
        },
        "workers": {
            wid: {
                "status": w.status.value,
                "current_task": w.current_task,
                "tasks_completed": w.tasks_completed,
                "container_id": w.container_id,  # NEW (may be None)
            }
            for wid, w in self._workers.items()
        },
        "levels": level_status["levels"],
        "is_complete": level_status["is_complete"],
    }
```

## Verification

```bash
python -c "
from zerg.orchestrator import Orchestrator
from zerg.config import ZergConfig
from zerg.launcher import ContainerLauncher

# Check orchestrator has updated methods
orch = Orchestrator('test', ZergConfig())

# Check status includes launcher_mode
import inspect
status_code = inspect.getsource(orch.status)
print('status() has launcher_mode:', 'launcher_mode' in status_code)

# Check _spawn_worker handles both modes
spawn_code = inspect.getsource(orch._spawn_worker)
print('_spawn_worker handles ContainerLauncher:', 'ContainerLauncher' in spawn_code or 'spawn_and_exec' in spawn_code)

# Check get_launcher_mode exists
print('get_launcher_mode exists:', hasattr(orch, 'get_launcher_mode'))
"
```

## Acceptance Criteria

- [ ] _spawn_worker() uses spawn_and_exec() for ContainerLauncher
- [ ] _spawn_worker() sets container_id on WorkerState
- [ ] _poll_workers() works for both launcher types
- [ ] _terminate_worker() uses launcher.terminate()
- [ ] status() includes launcher_mode
- [ ] status() includes container_id for workers
- [ ] Events include mode information
- [ ] No ruff errors: `ruff check zerg/orchestrator.py`
