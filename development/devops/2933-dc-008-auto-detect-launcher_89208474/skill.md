# DC-008: Add Auto-Detect Launcher Mode in Orchestrator

**Level**: 3 | **Critical Path**: No | **Estimate**: 20 min
**Dependencies**: DC-006

## Objective

Update `_create_launcher()` in `zerg/orchestrator.py` to auto-detect the appropriate launcher:
- Use ContainerLauncher if `.devcontainer/devcontainer.json` exists AND Docker image is built
- Otherwise use SubprocessLauncher
- Config override takes precedence

## Files Owned

- `zerg/orchestrator.py` (modify)

## Files to Read

- `zerg/launcher.py` (ContainerLauncher.docker_available, image_exists)
- `zerg/config.py` (WorkersConfig.launcher_type)

## Implementation

Update `_create_launcher()` method:

```python
from zerg.launcher import ContainerLauncher, LauncherConfig, LauncherType, SubprocessLauncher, WorkerLauncher


def _create_launcher(self) -> WorkerLauncher:
    """Create worker launcher based on config and environment.

    Priority:
    1. Explicit config override (subprocess or container)
    2. Auto-detect based on devcontainer presence and image availability

    Returns:
        Configured WorkerLauncher instance
    """
    # Get configured launcher type
    configured_type = self.config.workers.launcher_type

    # Build launcher config
    launcher_config = LauncherConfig(
        timeout_seconds=self.config.workers.timeout_minutes * 60,
        log_dir=Path(self.config.logging.directory),
    )

    # Check for explicit override
    if configured_type == "subprocess":
        logger.debug("Using SubprocessLauncher (explicit config)")
        launcher_config.launcher_type = LauncherType.SUBPROCESS
        return SubprocessLauncher(launcher_config)

    if configured_type == "container":
        logger.debug("Using ContainerLauncher (explicit config)")
        launcher_config.launcher_type = LauncherType.CONTAINER
        return ContainerLauncher(launcher_config)

    # Auto-detect mode
    return self._auto_detect_launcher(launcher_config)


def _auto_detect_launcher(self, config: LauncherConfig) -> WorkerLauncher:
    """Auto-detect appropriate launcher based on environment.

    Args:
        config: Launcher configuration

    Returns:
        SubprocessLauncher or ContainerLauncher
    """
    # Check for devcontainer
    devcontainer_path = self.repo_path / ".devcontainer" / "devcontainer.json"

    if not devcontainer_path.exists():
        logger.debug("No devcontainer found, using SubprocessLauncher")
        config.launcher_type = LauncherType.SUBPROCESS
        return SubprocessLauncher(config)

    # Check Docker availability
    if not ContainerLauncher.docker_available():
        logger.warning("Devcontainer exists but Docker not available, using SubprocessLauncher")
        config.launcher_type = LauncherType.SUBPROCESS
        return SubprocessLauncher(config)

    # Check if image is built
    if not ContainerLauncher.image_exists("zerg-worker"):
        logger.info("Devcontainer exists but image not built, using SubprocessLauncher")
        logger.info("Run 'zerg init --with-containers' to build the image")
        config.launcher_type = LauncherType.SUBPROCESS
        return SubprocessLauncher(config)

    # All conditions met for container mode
    logger.info("Using ContainerLauncher (auto-detected)")
    config.launcher_type = LauncherType.CONTAINER
    return ContainerLauncher(config)


def get_launcher_mode(self) -> str:
    """Get the current launcher mode.

    Returns:
        'subprocess', 'container', or 'auto'
    """
    if isinstance(self.launcher, ContainerLauncher):
        return "container"
    return "subprocess"
```

Also add a method to check if container mode is available:

```python
def container_mode_available(self) -> tuple[bool, str]:
    """Check if container mode can be used.

    Returns:
        Tuple of (available, reason)
    """
    devcontainer_path = self.repo_path / ".devcontainer" / "devcontainer.json"

    if not devcontainer_path.exists():
        return False, "No .devcontainer/devcontainer.json found"

    if not ContainerLauncher.docker_available():
        return False, "Docker is not available or not running"

    if not ContainerLauncher.image_exists("zerg-worker"):
        return False, "Docker image 'zerg-worker' not built (run zerg init --with-containers)"

    return True, "Container mode available"
```

## Verification

```bash
python -c "
from zerg.orchestrator import Orchestrator
from zerg.config import ZergConfig
from zerg.launcher import SubprocessLauncher, ContainerLauncher

# Test with default config (should auto-detect)
config = ZergConfig()
orch = Orchestrator('test-feature', config)

print(f'Launcher type: {type(orch.launcher).__name__}')
print(f'Mode: {orch.get_launcher_mode()}')

available, reason = orch.container_mode_available()
print(f'Container mode available: {available}')
print(f'Reason: {reason}')
"
```

## Acceptance Criteria

- [ ] _create_launcher() checks explicit config first
- [ ] Auto-detect checks devcontainer.json existence
- [ ] Auto-detect checks Docker availability
- [ ] Auto-detect checks image existence
- [ ] Falls back to subprocess if any condition fails
- [ ] get_launcher_mode() returns current mode
- [ ] container_mode_available() returns availability + reason
- [ ] Proper logging for mode selection
- [ ] No ruff errors: `ruff check zerg/orchestrator.py`
