# Shared Runtime State for Parallel Worktree Execution

## Problem

flowctl stores task runtime state (status, assignee, claimed_at, updated_at, evidence) in tracked files (`.flow/tasks/*.json`). This causes state fragmentation when using git worktrees for parallel execution:

1. `flowctl start` modifies `.flow/tasks/fn-1.1.json` in working directory
2. Worktree is created from branch → gets OLD state (status: todo)
3. Each worktree has its own copy of `.flow/tasks/`
4. Changes in one worktree don't propagate to others
5. `flowctl done` fails because task appears as `todo` in the execution context

This blocks parallel worker orchestration patterns where multiple agents work on independent tasks in isolated worktrees.

## Root Cause

Using a **tracked directory** (`.flow/tasks/*.json`) as **mutable runtime state** guarantees fragmentation — `git worktree add` and `git reset --hard` intentionally discard uncommitted tracked changes.

## Solution

Separate **definition** (tracked) from **runtime state** (shared, untracked):

| Data | Type | Location |
|------|------|----------|
| Task specs (`.md`) | Definition | `.flow/tasks/*.md` (tracked) |
| Task metadata (title, depends_on) | Definition | `.flow/tasks/*.json` (tracked) |
| Status, assignee, claimed_at, evidence | Runtime | `<state-dir>/tasks/*.state.json` (untracked) |

### State Directory Resolution

```python
def get_state_dir() -> Path:
    # 1. Explicit override (orchestrators can set)
    if state_dir := os.environ.get("FLOW_STATE_DIR"):
        return Path(state_dir).resolve()

    # 2. Git common-dir (shared across all worktrees automatically)
    common = subprocess.check_output(
        ["git", "rev-parse", "--git-common-dir", "--path-format=absolute"],
        text=True
    ).strip()
    return Path(common) / "flow-state"

    # 3. Fallback for non-git repos
    return get_flow_dir() / "state"
```

**Key insight:** `git rev-parse --git-common-dir` returns the shared `.git/` directory for all worktrees. This makes state automatically shared without env vars.

### StateStore Interface (Future-Proof)

```python
class StateStore(ABC):
    @abstractmethod
    def load_runtime(self, task_id: str) -> dict | None: ...

    @abstractmethod
    def save_runtime(self, task_id: str, data: dict) -> None: ...

    @abstractmethod
    def lock_task(self, task_id: str) -> ContextManager: ...

# Ship now:
class LocalFileStateStore(StateStore):
    """File-based state with fcntl locking."""
    pass

# Future (K8s, distributed):
# class HttpStateStore(StateStore): ...
# class RedisStateStore(StateStore): ...
```

### Runtime Fields

```python
RUNTIME_FIELDS = {
    "status",
    "updated_at",
    "claimed_at",
    "assignee",
    "claim_note",
    "evidence",
    "blocked_reason",
}
```

### Read Path (Merged View)

```python
def load_task(task_id: str) -> dict:
    # Load definition (tracked)
    def_path = get_flow_dir() / TASKS_DIR / f"{task_id}.json"
    definition = load_json(def_path)

    # Load runtime (state-dir)
    runtime = state_store.load_runtime(task_id)

    # Backward compat: if no runtime file, check for legacy fields in definition
    if runtime is None:
        runtime = {k: definition[k] for k in RUNTIME_FIELDS if k in definition}
        if not runtime:
            runtime = {"status": "todo"}

    # Merge: runtime overwrites definition for runtime fields
    return {**definition, **runtime}
```

### Write Path (Separate)

```python
def save_task_runtime(task_id: str, updates: dict) -> None:
    """Write runtime state only. Never touch definition file."""
    with state_store.lock_task(task_id):
        current = state_store.load_runtime(task_id) or {"status": "todo"}
        merged = {**current, **updates, "updated_at": now_iso()}
        state_store.save_runtime(task_id, merged)

# Definition writes remain separate (only for task create/edit commands)
def save_task_definition(task_id: str, definition: dict) -> None:
    """Write definition to tracked file."""
    def_path = get_flow_dir() / TASKS_DIR / f"{task_id}.json"
    # Filter out runtime fields
    clean_def = {k: v for k, v in definition.items() if k not in RUNTIME_FIELDS}
    atomic_write_json(def_path, clean_def)
```

## Backward Compatibility

| Scenario | Behavior |
|----------|----------|
| No state file exists | Read legacy runtime fields from definition |
| State file exists | Use state file, ignore legacy definition fields |
| Write operations | Always write to state-dir, never edit definition |
| Migration | Optional `flowctl migrate-state` command |

**No breaking changes for single-user flow-next.** Users who don't use worktrees see no difference.

## Concurrency

Parallel workers may update state simultaneously. Use file locking:

```python
# Per-task lock file
lock_path = state_dir / "locks" / f"{task_id}.lock"

with open(lock_path, 'w') as f:
    fcntl.flock(f, fcntl.LOCK_EX)
    try:
        # read-modify-write
    finally:
        fcntl.flock(f, fcntl.LOCK_UN)
```

## Commands Affected

| Command | Change |
|---------|--------|
| `flowctl start` | Write to state-dir |
| `flowctl done` | Write to state-dir |
| `flowctl block` | Write to state-dir |
| `flowctl show` | Merge definition + runtime |
| `flowctl tasks` | Merge definition + runtime |
| `flowctl ready` | Read merged state for dependency check |

## Debugging

Add helper command:

```bash
flowctl state-path
# Output: /path/to/.git/flow-state

flowctl state-path --task fn-1.1
# Output: /path/to/.git/flow-state/tasks/fn-1.1.state.json
```

## Future: Distributed Backends

The StateStore interface allows swapping backends via config:

```json
{
  "stateBackend": "local",
  "stateConfig": {
    "path": "/custom/path"
  }
}
```

Future backends (K8s, multi-machine):
- `http` — API server wrapping state
- `redis` — Shared Redis instance
- `s3` — For serverless/ephemeral workers

## Quick commands

```bash
# Run tests
python -m pytest plugins/flow-next/scripts/test_flowctl.py -v

# Manual verification
cd /tmp/test-repo && git init && .flow/bin/flowctl init
.flow/bin/flowctl state-path  # should show git-common-dir path
```

## Acceptance

- [ ] `flowctl start/done/block` write to state-dir, not definition file
- [ ] `flowctl show/tasks/ready` return merged definition + runtime
- [ ] State is shared across git worktrees (same git-common-dir)
- [ ] Existing repos without state-dir work (backward compat)
- [ ] Parallel `flowctl start` calls don't corrupt state (locking)
- [ ] `FLOW_STATE_DIR` env var overrides default
- [ ] `flowctl state-path` shows resolved path

## Tasks

1. **Implement StateStore interface and LocalFileStateStore** — Abstract interface, file-based implementation with locking
2. **Add state directory resolution** — git-common-dir default, env var override
3. **Split load_task into definition + runtime merge** — Backward-compatible read path
4. **Update write commands to use state-dir** — start, done, block
5. **Update read commands to merge state** — show, tasks, ready, list
6. **Add state-path debug command** — For troubleshooting
7. **Add migrate-state command** — Optional cleanup of legacy fields
8. **Update tests** — Mock state-dir, test worktree scenarios

## Non-Goals

- Changing task definition format
- Changing spec file format
- Auto-migrating existing repos (explicit opt-in)
- Distributed backends (future work)

## References

- Git worktree common-dir: `git rev-parse --git-common-dir`
- fcntl file locking: POSIX advisory locks
