# DC-006: Implement ContainerLauncher Base Class

**Level**: 1 | **Critical Path**: No | **Estimate**: 45 min

## Objective

Add `ContainerLauncher` class to `zerg/launcher.py` that implements the `WorkerLauncher` interface using Docker CLI commands.

## Files Owned

- `zerg/launcher.py` (modify)

## Implementation

Add after `SubprocessLauncher` class:

```python
class ContainerLauncher(WorkerLauncher):
    """Launch workers as Docker containers.

    Uses docker CLI to spawn isolated container environments for workers.
    Suitable for production deployments with resource isolation.
    """

    def __init__(self, config: LauncherConfig | None = None) -> None:
        """Initialize container launcher.

        Args:
            config: Launcher configuration
        """
        super().__init__(config)
        self._container_ids: dict[int, str] = {}
        self._image_name = "zerg-worker"

    def spawn(
        self,
        worker_id: int,
        feature: str,
        worktree_path: Path,
        branch: str,
        env: dict[str, str] | None = None,
    ) -> SpawnResult:
        """Spawn a new worker container.

        Args:
            worker_id: Unique worker identifier
            feature: Feature name being worked on
            worktree_path: Path to worker's git worktree
            branch: Git branch for worker
            env: Additional environment variables

        Returns:
            SpawnResult with handle or error
        """
        try:
            # Validate inputs
            if not isinstance(worker_id, int) or worker_id < 0:
                raise ValueError(f"Invalid worker_id: {worker_id}")

            container_name = f"zerg-worker-{worker_id}"

            # Build environment
            container_env = {
                "ZERG_WORKER_ID": str(worker_id),
                "ZERG_FEATURE": feature,
                "ZERG_BRANCH": branch,
            }

            # Add API key if available
            api_key = os.environ.get("ANTHROPIC_API_KEY")
            if api_key:
                container_env["ANTHROPIC_API_KEY"] = api_key

            # Validate and add additional env vars
            if self.config.env_vars:
                container_env.update(validate_env_vars(self.config.env_vars))
            if env:
                container_env.update(validate_env_vars(env))

            # Build docker run command
            cmd = [
                "docker", "run", "-d",
                "--name", container_name,
                "-v", f"{worktree_path.resolve()}:/workspace",
                "-w", "/workspace",
            ]

            # Add environment variables
            for key, value in container_env.items():
                cmd.extend(["-e", f"{key}={value}"])

            # Add network
            cmd.extend(["--network", "zerg-internal"])

            # Add image and command (keep alive)
            cmd.extend([self._image_name, "tail", "-f", "/dev/null"])

            # Start container
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=60,
            )

            if result.returncode != 0:
                return SpawnResult(
                    success=False,
                    worker_id=worker_id,
                    error=f"Docker run failed: {result.stderr}",
                )

            container_id = result.stdout.strip()[:12]
            self._container_ids[worker_id] = container_id

            # Create handle
            handle = WorkerHandle(
                worker_id=worker_id,
                container_id=container_id,
                status=WorkerStatus.INITIALIZING,
            )

            self._workers[worker_id] = handle
            logger.info(f"Spawned container {container_name} ({container_id})")

            return SpawnResult(success=True, worker_id=worker_id, handle=handle)

        except subprocess.TimeoutExpired:
            return SpawnResult(
                success=False,
                worker_id=worker_id,
                error="Docker run timed out",
            )
        except FileNotFoundError:
            return SpawnResult(
                success=False,
                worker_id=worker_id,
                error="Docker not found. Is Docker installed and running?",
            )
        except Exception as e:
            logger.error(f"Failed to spawn container for worker {worker_id}: {e}")
            return SpawnResult(success=False, worker_id=worker_id, error=str(e))

    def monitor(self, worker_id: int) -> WorkerStatus:
        """Check container status.

        Args:
            worker_id: Worker to check

        Returns:
            Current worker status
        """
        handle = self._workers.get(worker_id)
        container_id = self._container_ids.get(worker_id)

        if not handle or not container_id:
            return WorkerStatus.STOPPED

        try:
            result = subprocess.run(
                ["docker", "inspect", "-f", "{{.State.Status}}", container_id],
                capture_output=True,
                text=True,
                timeout=10,
            )

            if result.returncode != 0:
                handle.status = WorkerStatus.STOPPED
                return WorkerStatus.STOPPED

            state = result.stdout.strip()

            if state == "running":
                if handle.status == WorkerStatus.INITIALIZING:
                    handle.status = WorkerStatus.RUNNING
                return handle.status
            elif state == "exited":
                # Check exit code
                exit_result = subprocess.run(
                    ["docker", "inspect", "-f", "{{.State.ExitCode}}", container_id],
                    capture_output=True,
                    text=True,
                    timeout=10,
                )
                exit_code = int(exit_result.stdout.strip()) if exit_result.returncode == 0 else -1
                handle.exit_code = exit_code

                if exit_code == 0:
                    handle.status = WorkerStatus.STOPPED
                elif exit_code == 2:
                    handle.status = WorkerStatus.CHECKPOINTING
                else:
                    handle.status = WorkerStatus.CRASHED

                return handle.status
            else:
                return handle.status

        except Exception as e:
            logger.warning(f"Failed to monitor container {container_id}: {e}")
            return handle.status

    def terminate(self, worker_id: int, force: bool = False) -> bool:
        """Terminate a worker container.

        Args:
            worker_id: Worker to terminate
            force: Force termination with SIGKILL

        Returns:
            True if termination succeeded
        """
        container_id = self._container_ids.get(worker_id)
        handle = self._workers.get(worker_id)

        if not container_id or not handle:
            return False

        try:
            # Stop container
            cmd = ["docker", "kill" if force else "stop", container_id]
            subprocess.run(cmd, capture_output=True, timeout=30)

            # Remove container
            subprocess.run(
                ["docker", "rm", "-f", container_id],
                capture_output=True,
                timeout=10,
            )

            handle.status = WorkerStatus.STOPPED
            logger.info(f"Terminated container {container_id}")

            # Cleanup
            del self._container_ids[worker_id]
            return True

        except Exception as e:
            logger.error(f"Failed to terminate container {container_id}: {e}")
            return False

    def get_output(self, worker_id: int, tail: int = 100) -> str:
        """Get container logs.

        Args:
            worker_id: Worker to get output from
            tail: Number of lines from end

        Returns:
            Output string
        """
        container_id = self._container_ids.get(worker_id)
        if not container_id:
            return ""

        try:
            result = subprocess.run(
                ["docker", "logs", "--tail", str(tail), container_id],
                capture_output=True,
                text=True,
                timeout=10,
            )
            return result.stdout + result.stderr
        except Exception as e:
            return f"Error getting logs: {e}"

    @staticmethod
    def docker_available() -> bool:
        """Check if Docker is available and running."""
        try:
            result = subprocess.run(
                ["docker", "info"],
                capture_output=True,
                timeout=10,
            )
            return result.returncode == 0
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return False

    @staticmethod
    def image_exists(image_name: str = "zerg-worker") -> bool:
        """Check if the worker image exists."""
        try:
            result = subprocess.run(
                ["docker", "image", "inspect", image_name],
                capture_output=True,
                timeout=10,
            )
            return result.returncode == 0
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return False
```

## Verification

```bash
python -c "
from zerg.launcher import ContainerLauncher, WorkerLauncher, LauncherType
print(f'Inherits: {ContainerLauncher.__bases__}')
print(f'Methods: spawn, monitor, terminate, get_output')
print(f'LauncherType.CONTAINER: {LauncherType.CONTAINER.value}')
cl = ContainerLauncher()
print(f'Docker available: {cl.docker_available()}')
"
```

## Acceptance Criteria

- [ ] ContainerLauncher inherits from WorkerLauncher
- [ ] spawn() creates container with docker run
- [ ] monitor() checks container state via docker inspect
- [ ] terminate() stops and removes container
- [ ] get_output() returns docker logs
- [ ] docker_available() static method works
- [ ] image_exists() static method works
- [ ] No ruff errors: `ruff check zerg/launcher.py`
