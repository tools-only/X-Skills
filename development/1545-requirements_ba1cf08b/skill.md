# Feature Requirements: container-execution

## Metadata
- **Feature**: container-execution
- **Status**: APPROVED
- **Approved**: 2026-01-28
- **Created**: 2026-01-28
- **Author**: ZERG Plan Mode

---

## 1. Problem Statement

### 1.1 Background
ZERG workers currently only run as local subprocesses. The ContainerLauncher class is fully implemented (~550 LOC) but has never been exercised end-to-end. The Docker image (`zerg-worker`) has never been built. The existing plan (`claudedocs/plan-container-dogfooding.md`) outlines build-and-verify steps but doesn't cover: resource limits, OAuth credential passthrough, health checks, cross-platform compatibility, or production hardening.

### 1.2 Problem
- Docker image doesn't exist — `docker images | grep zerg-worker` returns nothing
- ContainerLauncher has never been run against real Docker — only unit-tested with mocked subprocess calls
- No resource limits on containers — a runaway worker could consume all host memory/CPU
- OAuth credential mounting is implemented but untested — users without API keys can't use container mode
- No health check or auto-recovery for container workers
- No log collection from container workers to host
- `docker-compose.yaml` defines resource limits (4G/2CPU) but ContainerLauncher bypasses compose entirely
- `privileged: true` in docker-compose.yaml is a security concern

### 1.3 Impact
Without container mode, ZERG cannot provide worker isolation. Subprocess workers share the host filesystem, environment, and can interfere with each other or the user's system. Container mode is the security boundary that makes parallel execution safe.

---

## 2. Users

### 2.1 Primary Users
Developers running ZERG on macOS (Docker Desktop, ARM64) or Linux (native Docker) to execute parallel feature builds with container isolation.

### 2.2 User Stories
- As a developer, I want to run `zerg rush --mode container` and have workers execute in isolated Docker containers so that parallel workers can't interfere with each other
- As a developer, I want container workers to authenticate with Claude Code via my existing OAuth session so I don't need to manage API keys
- As a developer, I want container workers to have resource limits so a runaway worker can't crash my machine
- As a developer, I want to see container worker logs via `zerg logs` so I can debug failures

---

## 3. Functional Requirements

### 3.1 Core Capabilities

| ID | Requirement | Priority | Notes |
|----|-------------|----------|-------|
| FR-001 | Docker image builds successfully on macOS ARM64 and Linux x86_64 | Must | Multi-platform support |
| FR-002 | `zerg rush --mode container --workers N` spawns N container workers | Must | Core execution path |
| FR-003 | Container workers execute tasks, commit changes, report status | Must | Full task lifecycle |
| FR-004 | Resource limits (memory, CPU) applied per container | Must | Default: 4G RAM, 2 CPU |
| FR-005 | ANTHROPIC_API_KEY passthrough to containers | Must | Primary auth method |
| FR-006 | OAuth credential passthrough via ~/.claude mount | Must | Secondary auth method |
| FR-007 | Git worktree operations work inside containers | Must | Commits, branches, pushes |
| FR-008 | Worker state (RUNNING/STOPPED/CRASHED) visible via `zerg status` | Must | Uses state file IPC (just implemented) |
| FR-009 | Container cleanup on worker exit (success, failure, or crash) | Must | No orphan containers |
| FR-010 | `zerg rush --mode auto` detects image and selects container mode | Should | Falls back to subprocess |
| FR-011 | Container worker logs accessible via `docker logs` and `zerg logs` | Should | Log collection |
| FR-012 | Health check detects stuck containers and marks them CRASHED | Should | Timeout-based |
| FR-013 | `zerg build` command to build the Docker image | Could | Convenience wrapper |
| FR-014 | Remove `privileged: true` from docker-compose.yaml | Must | Security hardening |

### 3.2 Inputs
- `zerg rush --mode container --workers N --feature <name>`
- ANTHROPIC_API_KEY or OAuth credentials in ~/.claude
- Docker daemon running and accessible

### 3.3 Outputs
- N running containers named `zerg-worker-0..N-1`
- Task results committed to worker branches in worktrees
- Worker state in `.zerg/state/{feature}.json`
- Container logs accessible via `docker logs zerg-worker-N`

### 3.4 Business Rules
- Container workers MUST NOT have access to ~/.ssh, ~/.aws, or ~/.config (only ~/.claude)
- Container workers MUST run as non-root (current UID/GID mapping already implemented)
- Containers MUST be cleaned up even on orchestrator crash (use `--rm` or explicit cleanup)

---

## 4. Non-Functional Requirements

### 4.1 Performance
- Container spawn: <30s per worker (image pull excluded)
- Container overhead vs subprocess: <5% task execution time penalty
- Resource limits: 4G RAM, 2 CPU per container (configurable via .zerg/config.yaml)

### 4.2 Security
- No privileged mode
- Non-root execution (UID/GID pass-through)
- Mount only: worktree, .zerg/state, ~/.claude, .git (for commits)
- No host network mode — use bridge network
- Environment variable allowlist enforced (existing validate_env_vars)

### 4.3 Reliability
- Orphan container detection and cleanup on orchestrator start
- Container health check via `docker inspect` in poll loop
- Graceful shutdown: SIGTERM → 10s wait → SIGKILL (already implemented)

### 4.4 Platform Support
- macOS Docker Desktop (ARM64) — primary
- Linux native Docker (x86_64) — secondary
- No Windows support required

---

## 5. Scope

### 5.1 In Scope
- Build and verify Docker image
- E2E dogfood: plan → design → rush --mode container → status → merge
- Resource limits on containers
- OAuth credential passthrough
- Cross-platform Dockerfile (multi-arch)
- Fix security issues (remove privileged, constrain mounts)
- Integration tests with real Docker (gated on CI Docker availability)
- Log collection from containers
- Health check and auto-recovery for stuck workers

### 5.2 Out of Scope
- Docker-in-Docker (Workflow B: container spawning containers) — deferred
- docker-compose orchestration — keep using docker run
- wait-for-services.sh (postgres/redis profiles) — separate feature
- Dynamic devcontainer generation from devcontainer_features.py — separate feature
- Kubernetes/remote Docker host support — future

### 5.3 Assumptions
- Docker daemon is running and accessible via `docker` CLI
- User has built the image (or `zerg build` will build it)
- ANTHROPIC_API_KEY is set OR user has valid OAuth session in ~/.claude

### 5.4 Constraints
- Must not break existing subprocess mode
- Must use existing ContainerLauncher API (no interface changes)
- Must work with existing state file IPC (just implemented)

---

## 6. Dependencies

### 6.1 Internal Dependencies
| Dependency | Type | Status |
|------------|------|--------|
| State file IPC | Required | Complete (a189fc7) |
| ContainerLauncher class | Required | Implemented, untested E2E |
| SubprocessLauncher | Required | Working, must not regress |
| worker_entry.sh | Required | Implemented, untested in real container |

### 6.2 External Dependencies
| Dependency | Type | Owner |
|------------|------|-------|
| Docker daemon | Required | User |
| Docker Desktop (macOS) | Required | User |
| ANTHROPIC_API_KEY or OAuth | Required | User |
| Claude Code CLI (npm) | Required | Dockerfile |

---

## 7. Acceptance Criteria

### 7.1 Definition of Done
- [ ] Docker image builds on macOS ARM64
- [ ] Docker image builds on Linux x86_64
- [ ] `zerg rush --mode container --workers 2 --dry-run` completes
- [ ] `zerg rush --mode container --workers 2` executes a real feature
- [ ] Workers write state visible via `zerg status` in another terminal
- [ ] Workers commit code to their branches
- [ ] Merge succeeds after level completion
- [ ] Resource limits applied (verify via `docker inspect`)
- [ ] No orphan containers after rush completes
- [ ] OAuth authentication works in containers
- [ ] All existing tests pass (no subprocess regression)
- [ ] New integration tests for container mode pass
- [ ] privileged: true removed from docker-compose.yaml

### 7.2 Test Scenarios

| ID | Scenario | Given | When | Then |
|----|----------|-------|------|------|
| TC-001 | Image build | Dockerfile exists | `docker build -t zerg-worker .` | Image built, claude --version works |
| TC-002 | Container spawn | Image exists | `zerg rush --mode container --dry-run` | Dry run shows container mode |
| TC-003 | Full E2E | Feature planned+designed | `zerg rush --mode container --workers 2` | Tasks complete, branches merged |
| TC-004 | State visibility | Workers running | `zerg status` from another terminal | Worker data appears |
| TC-005 | Resource limits | Container running | `docker inspect zerg-worker-0` | Memory=4G, CPU=2 |
| TC-006 | Cleanup on crash | Worker crashes | Orchestrator poll loop runs | Container removed, state=CRASHED |
| TC-007 | OAuth auth | ~/.claude has valid session | Container worker starts | Claude Code authenticates |
| TC-008 | API key auth | ANTHROPIC_API_KEY set | Container worker starts | Claude Code authenticates |
| TC-009 | No orphans | Rush completes | `docker ps -a \| grep zerg-worker` | No containers left |
| TC-010 | Subprocess unaffected | No Docker | `zerg rush --mode subprocess` | Works as before |

### 7.3 Success Metrics
- Container workers complete tasks at same success rate as subprocess workers
- Zero orphan containers after any rush completion
- `zerg status` shows live worker data during container execution

---

## 8. Open Questions

| ID | Question | Status |
|----|----------|--------|
| Q-001 | Should `zerg build` auto-detect platform and build multi-arch? | Open |
| Q-002 | Should resource limits be configurable in .zerg/config.yaml? | Open (default to yes) |
| Q-003 | Should we add `--auto-build` flag to `zerg rush` to build image if missing? | Open |

---

## 9. Approval

| Role | Name | Date | Signature |
|------|------|------|-----------|
| Product | | | PENDING |
| Engineering | | | PENDING |
