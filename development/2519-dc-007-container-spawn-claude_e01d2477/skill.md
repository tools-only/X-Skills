# DC-007: Add Container Spawn with Claude Execution

**Level**: 2 | **Critical Path**: No | **Estimate**: 40 min
**Dependencies**: DC-006

## Objective

Extend `ContainerLauncher` to execute Claude inside containers via a worker entry script. Create `.zerg/worker_entry.sh` that invokes `claude` with proper flags.

## Files Owned

- `zerg/launcher.py` (modify - add exec_claude method)
- `.zerg/worker_entry.sh` (create)

## Implementation

### 1. Create `.zerg/worker_entry.sh`

```bash
#!/bin/bash
# ZERG Worker Entry Script
# Invokes Claude Code CLI inside container

set -e

WORKER_ID=${ZERG_WORKER_ID:-0}
FEATURE=${ZERG_FEATURE:-"unknown"}
BRANCH=${ZERG_BRANCH:-"main"}

echo "================================================"
echo "ZERG Worker $WORKER_ID starting..."
echo "Feature: $FEATURE"
echo "Branch: $BRANCH"
echo "Workspace: $(pwd)"
echo "================================================"

# Check Claude is available
if ! command -v claude &> /dev/null; then
    echo "ERROR: Claude Code CLI not found"
    echo "Installing Claude Code..."
    npm install -g @anthropic-ai/claude-code
fi

# Verify API key
if [ -z "$ANTHROPIC_API_KEY" ]; then
    echo "ERROR: ANTHROPIC_API_KEY not set"
    exit 1
fi

# State file for communication with orchestrator
STATE_FILE="/workspace/.zerg/state/worker-${WORKER_ID}.json"
mkdir -p "$(dirname "$STATE_FILE")"

# Write initial state
echo "{\"worker_id\": $WORKER_ID, \"status\": \"running\", \"started_at\": \"$(date -Iseconds)\"}" > "$STATE_FILE"

# Launch Claude Code
# --dangerously-skip-permissions: Required for non-interactive execution
# STDIN is provided via task file or prompt
echo "Launching Claude Code..."

claude --dangerously-skip-permissions \
    --output-format json \
    2>&1 | tee "/workspace/.zerg/logs/worker-${WORKER_ID}.log"

EXIT_CODE=$?

# Update state with exit
echo "{\"worker_id\": $WORKER_ID, \"status\": \"exited\", \"exit_code\": $EXIT_CODE, \"ended_at\": \"$(date -Iseconds)\"}" > "$STATE_FILE"

echo "Worker $WORKER_ID exited with code $EXIT_CODE"
exit $EXIT_CODE
```

### 2. Update ContainerLauncher in `zerg/launcher.py`

Add method to execute claude in container:

```python
class ContainerLauncher(WorkerLauncher):
    # ... existing code ...

    def exec_claude(
        self,
        worker_id: int,
        prompt: str | None = None,
        prompt_file: Path | None = None,
    ) -> bool:
        """Execute Claude inside a running container.

        Args:
            worker_id: Worker container to exec in
            prompt: Direct prompt string
            prompt_file: Path to prompt file (relative to workspace)

        Returns:
            True if exec started successfully
        """
        container_id = self._container_ids.get(worker_id)
        if not container_id:
            logger.error(f"No container for worker {worker_id}")
            return False

        try:
            # Execute worker entry script
            cmd = [
                "docker", "exec", "-d",
                container_id,
                "/workspace/.zerg/worker_entry.sh",
            ]

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=30,
            )

            if result.returncode != 0:
                logger.error(f"Failed to exec in container: {result.stderr}")
                return False

            # Update handle status
            handle = self._workers.get(worker_id)
            if handle:
                handle.status = WorkerStatus.RUNNING

            logger.info(f"Started Claude in container {container_id}")
            return True

        except Exception as e:
            logger.error(f"Failed to exec claude in container: {e}")
            return False

    def spawn_and_exec(
        self,
        worker_id: int,
        feature: str,
        worktree_path: Path,
        branch: str,
        env: dict[str, str] | None = None,
    ) -> SpawnResult:
        """Spawn container and immediately start Claude.

        Convenience method that combines spawn() + exec_claude().

        Args:
            worker_id: Unique worker identifier
            feature: Feature name
            worktree_path: Path to worktree
            branch: Git branch
            env: Additional environment

        Returns:
            SpawnResult
        """
        # First spawn the container
        result = self.spawn(worker_id, feature, worktree_path, branch, env)

        if not result.success:
            return result

        # Wait for container to be ready
        import time
        for _ in range(10):
            status = self.monitor(worker_id)
            if status == WorkerStatus.RUNNING:
                break
            time.sleep(0.5)

        # Execute Claude
        if not self.exec_claude(worker_id):
            return SpawnResult(
                success=False,
                worker_id=worker_id,
                error="Container started but failed to exec Claude",
            )

        return result

    def wait_for_completion(
        self,
        worker_id: int,
        timeout: float = 3600,
        poll_interval: float = 5.0,
    ) -> WorkerStatus:
        """Wait for worker to complete.

        Args:
            worker_id: Worker to wait for
            timeout: Maximum wait time in seconds
            poll_interval: Time between status checks

        Returns:
            Final worker status
        """
        import time
        start = time.time()

        while time.time() - start < timeout:
            status = self.monitor(worker_id)
            if status in (WorkerStatus.STOPPED, WorkerStatus.CRASHED, WorkerStatus.CHECKPOINTING):
                return status
            time.sleep(poll_interval)

        return WorkerStatus.RUNNING  # Still running after timeout
```

## Verification

```bash
# Check script exists and is executable
test -x .zerg/worker_entry.sh && echo "Script executable: OK"

# Check script contains claude invocation
grep -q 'claude --dangerously-skip-permissions' .zerg/worker_entry.sh && echo "Claude invocation: OK"

# Check ContainerLauncher has new methods
python -c "
from zerg.launcher import ContainerLauncher
cl = ContainerLauncher()
print('exec_claude:', hasattr(cl, 'exec_claude'))
print('spawn_and_exec:', hasattr(cl, 'spawn_and_exec'))
print('wait_for_completion:', hasattr(cl, 'wait_for_completion'))
"
```

## Acceptance Criteria

- [ ] `.zerg/worker_entry.sh` exists and is executable (chmod +x)
- [ ] Script checks for Claude CLI, API key
- [ ] Script writes state to `/workspace/.zerg/state/worker-{id}.json`
- [ ] ContainerLauncher.exec_claude() runs script in container
- [ ] ContainerLauncher.spawn_and_exec() combines spawn + exec
- [ ] wait_for_completion() polls until done
- [ ] No ruff errors
