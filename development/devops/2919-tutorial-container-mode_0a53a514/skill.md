# Tutorial: Container Mode

This tutorial covers running ZERG workers inside Docker containers. Container mode provides isolated, reproducible execution environments where each worker runs in its own container with controlled resources, network access, and filesystem mounts.

Use container mode when you want:

- Isolation between workers (no shared process space)
- Reproducible environments across machines
- Resource limits enforced by Docker (CPU, memory)
- A clean separation between worker filesystems and the host

## Table of Contents

- [Prerequisites](#prerequisites)
- [Building the Worker Image](#building-the-worker-image)
- [Configuration](#configuration)
- [Authentication](#authentication)
- [Launching in Container Mode](#launching-in-container-mode)
- [Volume Mounts and State Sharing](#volume-mounts-and-state-sharing)
- [Resource Limits](#resource-limits)
- [Debugging Containers](#debugging-containers)
- [Troubleshooting](#troubleshooting)

---

## Prerequisites

Before using container mode, confirm the following:

```bash
# Docker must be installed and running
docker --version
# Expected: Docker version 24.x or later

# Verify Docker daemon is accessible
docker info > /dev/null 2>&1 && echo "Docker is running" || echo "Docker is NOT running"

# Confirm you can run containers
docker run --rm hello-world
```

You also need a completed design phase. Container mode is a launch option for `/zerg:rush`, not a separate workflow. Complete `/zerg:plan` and `/zerg:design` first. See [[Tutorial-Minerals-Store]] for a full walkthrough of those phases.

---

## Building the Worker Image

ZERG workers run inside a Docker image named `zerg-worker` by default. Build it from the project root:

```bash
docker build -t zerg-worker -f .devcontainer/Dockerfile .
```

The image should include:

- Python 3.12 or later
- Node.js (for Claude Code)
- Git
- The `claude` CLI tool
- Any project-specific dependencies

Verify the image exists:

```bash
docker images zerg-worker
```

Expected output:

```
REPOSITORY    TAG       IMAGE ID       CREATED        SIZE
zerg-worker   latest    a1b2c3d4e5f6   2 hours ago    1.2GB
```

If you are using a custom image name, pass it to the launcher via configuration (see the Configuration section below).

---

## Configuration

Container mode is controlled through `.zerg/config.yaml`. The relevant settings are under the `workers` and `resources` sections:

```yaml
# .zerg/config.yaml

workers:
  max_concurrent: 5
  timeout_minutes: 60
  retry_attempts: 2
  launcher_type: container    # Set this to "container" instead of "subprocess"

resources:
  cpu_cores: 2
  memory_gb: 4
  container_memory_limit: "4g"
  container_cpu_limit: 2.0

security:
  container_readonly: true    # Mount root filesystem as read-only
```

Key settings:

| Setting | Default | Description |
|---------|---------|-------------|
| `workers.launcher_type` | `subprocess` | Set to `container` for Docker execution |
| `resources.container_memory_limit` | `4g` | Docker `--memory` flag value |
| `resources.container_cpu_limit` | `2.0` | Docker `--cpus` flag value |
| `security.container_readonly` | `true` | Mount container root as read-only |

You can also activate container mode per-run without changing config:

```
/zerg:rush --workers=3 --mode container
```

---

## Authentication

Container workers need to authenticate with the Claude API. ZERG supports two methods.

### Method 1: OAuth (Claude Pro/Team Accounts)

OAuth authentication works by mounting your local `~/.claude` directory into the container. This is the default method when the directory exists on the host.

ZERG handles this automatically. When `ContainerLauncher` detects `~/.claude` on the host, it adds these volume mounts:

```
-v ~/.claude:/home/worker/.claude
-v ~/.claude.json:/home/worker/.claude.json
-e HOME=/home/worker
```

**No additional configuration is required.** If you are logged into Claude Code on your host machine, container workers inherit that session.

To verify your OAuth credentials are present:

```bash
ls -la ~/.claude/
# Should contain session files

ls -la ~/.claude.json
# Should exist and contain OAuth tokens
```

**Security note**: The `~/.claude` mount is read-write because Claude Code writes debug logs and session state to this directory. The mount shares your host credentials with all worker containers. This is appropriate for local development. Do not use OAuth mounts in shared or CI environments where other users could access the containers.

### Method 2: API Key

API key authentication passes your `ANTHROPIC_API_KEY` as an environment variable into each container. Use this method in CI/CD pipelines or when OAuth is not available.

Set the API key in your environment before launching:

```bash
export ANTHROPIC_API_KEY="sk-ant-api03-your-key-here"
```

Or place it in a `.env` file in the project root:

```
ANTHROPIC_API_KEY=sk-ant-api03-your-key-here
```

ZERG reads the key from either source. The `ContainerLauncher` checks `os.environ` first, then falls back to parsing `.env`:

```python
# From zerg/launcher.py (simplified)
api_key = os.environ.get("ANTHROPIC_API_KEY")
if not api_key:
    env_file = Path(".env")
    if env_file.exists():
        for line in env_file.read_text().splitlines():
            if line.startswith("ANTHROPIC_API_KEY"):
                _, _, val = line.partition("=")
                api_key = val.strip().strip("'\"")
```

The key is passed as `-e ANTHROPIC_API_KEY=...` to `docker run`. It is not baked into the image and does not persist in image layers.

**Security note**: Treat `.env` files containing API keys as sensitive. Add `.env` to your `.gitignore` and `.dockerignore` to prevent accidental commits or image inclusion.

### Choosing a Method

| Scenario | Recommended Method |
|----------|--------------------|
| Local development, logged into Claude Code | OAuth |
| CI/CD pipeline | API Key |
| Shared development machine | API Key (avoid sharing OAuth tokens) |
| Air-gapped environment | API Key (no OAuth server access) |

---

## Launching in Container Mode

With configuration and authentication in place, launch workers:

```
/zerg:rush --workers=3 --mode container
```

ZERG performs these steps:

1. Builds git worktrees for each worker (same as subprocess mode)
2. Removes any existing containers with matching names
3. Starts a Docker container per worker with volume mounts and environment variables
4. Waits for each container to become ready (up to 30 seconds)
5. Verifies the worker process started inside the container (up to 120 seconds)
6. Starts the orchestrator to coordinate level transitions

Container output during launch:

```
Container mode: using image zerg-worker
Launching Worker 0: zerg-worker-0 (memory=4g, cpus=2.0)
Launching Worker 1: zerg-worker-1 (memory=4g, cpus=2.0)
Launching Worker 2: zerg-worker-2 (memory=4g, cpus=2.0)
All containers running. Orchestrator started.
```

Verify containers are running:

```bash
docker ps --filter "name=zerg-worker" --format "table {{.Names}}\t{{.Status}}\t{{.Ports}}"
```

```
NAMES            STATUS          PORTS
zerg-worker-0    Up 2 minutes
zerg-worker-1    Up 2 minutes
zerg-worker-2    Up 2 minutes
```

---

## Volume Mounts and State Sharing

Each container receives several volume mounts to function correctly. Understanding these mounts is important for debugging.

### Workspace Mount

```
-v /path/to/.zerg-worktrees/feature/worker-0:/workspace
```

The worker's git worktree is mounted at `/workspace` inside the container. This is the container's working directory. All file creation and modification happens here.

### State Directory Mount

```
-v /path/to/repo/.zerg/state:/workspace/.zerg/state
```

The shared state directory is mounted from the main repository. This allows the orchestrator (running on the host) and workers (running in containers) to read and write the same state files. The state file tracks task status, worker assignments, and level completion.

### Git Mounts

Git worktrees require access to the main repository's `.git` directory for object storage and refs. ZERG mounts two paths:

```
-v /path/to/repo/.git:/repo/.git
-v /path/to/repo/.git/worktrees/worker-0:/workspace/.git-worktree
```

The worker entry script (`worker_entry.sh`) reconfigures git paths inside the container so that `git commit` and `git push` work against the correct repository.

Environment variables passed to support this:

```
ZERG_GIT_WORKTREE_DIR=/workspace/.git-worktree
ZERG_GIT_MAIN_DIR=/repo/.git
```

### OAuth Mounts (if present)

```
-v ~/.claude:/home/worker/.claude
-v ~/.claude.json:/home/worker/.claude.json
```

See the Authentication section above.

### Mount Summary

| Host Path | Container Path | Purpose |
|-----------|---------------|---------|
| `.zerg-worktrees/{feature}/worker-N` | `/workspace` | Worker's git worktree |
| `{repo}/.zerg/state` | `/workspace/.zerg/state` | Shared task state |
| `{repo}/.git` | `/repo/.git` | Git object store |
| `{repo}/.git/worktrees/worker-N` | `/workspace/.git-worktree` | Worktree metadata |
| `~/.claude` | `/home/worker/.claude` | OAuth credentials |
| `~/.claude.json` | `/home/worker/.claude.json` | OAuth tokens |

---

## Resource Limits

Container mode enforces resource limits through Docker flags. This prevents a single worker from consuming all host resources.

### Memory

Set via `resources.container_memory_limit` in config or the `--memory` Docker flag:

```yaml
resources:
  container_memory_limit: "4g"
```

If a worker exceeds this limit, Docker kills the container with an OOM (out of memory) error. The orchestrator detects the crash and can retry the task.

Recommended values:

| Workload | Memory Limit |
|----------|-------------|
| Small tasks (types, configs) | 2g |
| Medium tasks (services, logic) | 4g |
| Large tasks (full test suites) | 8g |

### CPU

Set via `resources.container_cpu_limit` in config or the `--cpus` Docker flag:

```yaml
resources:
  container_cpu_limit: 2.0
```

A value of `2.0` means the container can use up to 2 CPU cores. Docker throttles the container if it exceeds this limit (no kill, just slowdown).

### Calculating Host Capacity

If you are running 5 workers, each with 4g memory and 2 CPUs:

```
Total memory needed: 5 * 4g = 20g
Total CPUs needed:   5 * 2  = 10 cores
```

Add overhead for the orchestrator and host OS. A machine with 32GB RAM and 12 cores comfortably runs 5 workers at these limits.

### Per-Container Override

To override limits for a specific run without changing config:

```bash
# In .zerg/config.yaml, set defaults
# Then override at launch time via environment
ZERG_CONTAINER_MEMORY=8g ZERG_CONTAINER_CPUS=4.0 /zerg:rush --workers=2 --mode container
```

---

## Debugging Containers

When a worker fails inside a container, you need to inspect the container state.

### View Container Logs

```bash
# Stream logs from a specific worker
docker logs zerg-worker-0

# Follow logs in real time
docker logs -f zerg-worker-0

# Show last 50 lines
docker logs --tail 50 zerg-worker-0
```

### Execute a Shell Inside a Running Container

```bash
docker exec -it zerg-worker-0 bash
```

Once inside, inspect the workspace:

```bash
# Check working directory
pwd
# Expected: /workspace

# Verify git is configured
git status
git log --oneline -5

# Check environment variables
env | grep ZERG

# Look at the spec files
cat .gsd/specs/minerals-store/task-graph.json | python3 -m json.tool

# Check if Claude Code is available
claude --version
```

### Inspect a Stopped Container

If the container has exited, you can still read its filesystem:

```bash
# Check exit code
docker inspect zerg-worker-0 --format '{{.State.ExitCode}}'

# Check OOM status
docker inspect zerg-worker-0 --format '{{.State.OOMKilled}}'

# Copy logs out of the container
docker cp zerg-worker-0:/workspace/.zerg/logs/. ./debug-logs/

# Start the stopped container for inspection
docker start zerg-worker-0
docker exec -it zerg-worker-0 bash
```

### Worker Exit Codes

The exit code tells you why the worker stopped:

| Exit Code | Meaning | Action |
|-----------|---------|--------|
| 0 | All tasks completed | No action needed |
| 1 | Unrecoverable error | Check `docker logs` for stack trace |
| 2 | Context limit reached (70%) | Orchestrator restarts automatically |
| 3 | All remaining tasks blocked | Check task dependencies and failures |
| 130 | Received stop signal | Normal shutdown via `/zerg:stop` |
| 137 | Killed by Docker (OOM) | Increase `container_memory_limit` |

### Reading Worker-Specific Logs

Worker logs are written to the shared state directory:

```bash
# Host-side log access
cat .zerg/logs/workers/worker-0.stdout.log

# Progress log (all workers write to this)
cat .gsd/specs/minerals-store/progress.md
```

### Container Resource Usage

Monitor real-time resource consumption:

```bash
# All ZERG containers
docker stats --filter "name=zerg-worker" --format "table {{.Name}}\t{{.CPUPerc}}\t{{.MemUsage}}\t{{.MemPerc}}"
```

```
NAME             CPU %   MEM USAGE / LIMIT   MEM %
zerg-worker-0    45.2%   1.8GiB / 4GiB       45.0%
zerg-worker-1    38.7%   2.1GiB / 4GiB       52.5%
zerg-worker-2    12.3%   0.9GiB / 4GiB       22.5%
```

---

## Troubleshooting

### "Failed to start container"

**Cause**: Docker daemon is not running, or the image does not exist.

```bash
# Check Docker is running
docker info

# Check image exists
docker images zerg-worker

# Rebuild if missing
docker build -t zerg-worker -f .devcontainer/Dockerfile .
```

### "Container failed to become ready"

**Cause**: The container started but did not reach a running state within 30 seconds. Often caused by a missing entrypoint or broken image.

```bash
# Check container status
docker ps -a --filter "name=zerg-worker-0"

# Check container logs for startup errors
docker logs zerg-worker-0
```

### "Worker process failed to start"

**Cause**: The container is running, but the worker entry script failed. This has a 120-second timeout.

```bash
# Check if the entry script exists in the worktree
ls -la ../.zerg-worktrees/minerals-store/worker-0/.zerg/worker_entry.sh

# Exec into the container and run the entry script manually
docker exec -it zerg-worker-0 bash
cat /workspace/.zerg/worker_entry.sh
bash /workspace/.zerg/worker_entry.sh
```

### OOM Killed (Exit Code 137)

**Cause**: The worker exceeded the container memory limit.

```bash
# Confirm OOM
docker inspect zerg-worker-0 --format '{{.State.OOMKilled}}'
# Expected: true
```

**Fix**: Increase the memory limit in `.zerg/config.yaml`:

```yaml
resources:
  container_memory_limit: "8g"
```

### Authentication Failures

**Symptom**: Worker logs show "Authentication failed" or "API key not found".

For OAuth:

```bash
# Verify ~/.claude exists and contains session data
ls -la ~/.claude/

# Check that the mount is present in the container
docker exec zerg-worker-0 ls -la /home/worker/.claude/
```

For API key:

```bash
# Verify the key is set
echo $ANTHROPIC_API_KEY | head -c 20

# Check it reached the container
docker exec zerg-worker-0 env | grep ANTHROPIC_API_KEY
```

### Git Errors Inside Container

**Symptom**: Worker logs show "fatal: not a git repository" or similar git errors.

```bash
# Check the git worktree configuration inside the container
docker exec zerg-worker-0 cat /workspace/.git
# Should point to the worktree metadata

# Verify the main .git is accessible
docker exec zerg-worker-0 ls /repo/.git/HEAD
```

If the `.git` file in `/workspace` has a stale path, the entry script should fix it. Check:

```bash
docker exec zerg-worker-0 env | grep ZERG_GIT
# Should show:
# ZERG_GIT_WORKTREE_DIR=/workspace/.git-worktree
# ZERG_GIT_MAIN_DIR=/repo/.git
```

### Cleaning Up Stale Containers

If a previous run left containers behind:

```bash
# Stop and remove all ZERG worker containers
docker rm -f $(docker ps -aq --filter "name=zerg-worker")

# Remove worktrees
rm -rf ../.zerg-worktrees/minerals-store/
```

Or use the built-in cleanup:

```
/zerg:cleanup minerals-store
```

---

## Summary

Container mode provides isolated execution for ZERG workers with enforced resource limits and reproducible environments. The key points:

- Set `launcher_type: container` in config or pass `--mode container` at launch
- Build the `zerg-worker` Docker image before first use
- Authentication works via OAuth (`~/.claude` mount) or API key (`ANTHROPIC_API_KEY` env var)
- Each container receives volume mounts for the workspace, git metadata, shared state, and OAuth credentials
- Resource limits (memory, CPU) are enforced by Docker and configurable per-project
- Debug failed containers with `docker logs`, `docker exec`, and `docker inspect`
- Worker exit codes (0, 1, 2, 3, 137) tell you exactly why a worker stopped

Previous tutorial: [[Tutorial-Minerals-Store]] for the full ZERG workflow walkthrough.
