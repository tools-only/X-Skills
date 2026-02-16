# Tuning Guide

This page explains when and how to adjust ZERG configuration values for different project sizes, hardware constraints, and reliability requirements. All options referenced here are documented in [[Configuration]].

---

## Worker Count

The `workers.max_concurrent` setting controls how many Claude Code instances run in parallel. This is the single most impactful tuning parameter.

### Choosing a Value

| Project Size | Files per Feature | Recommended Workers | Reasoning |
|-------------|-------------------|--------------------:|-----------|
| Small | < 20 files | 2-3 | Overhead of coordination exceeds parallelism gains |
| Medium | 20-50 files | 4-6 | Good balance of parallelism and resource usage |
| Large | 50-100 files | 6-10 | Significant speedup if tasks are well-partitioned |
| Very Large | > 100 files | 8-15 | Diminishing returns above 10-12 in practice |

### Constraints to Consider

**API rate limits.** Each worker makes independent API calls. If you are on a plan with rate limits, reduce workers to stay within your quota.

**System memory.** Each subprocess worker consumes approximately 200-400 MB of memory. Container workers use whatever `resources.container_memory_limit` specifies. A machine with 16 GB of RAM can comfortably run 5-8 subprocess workers alongside normal development tools.

**CPU cores.** Workers are I/O-bound (waiting on API responses) more than CPU-bound. You can safely run more workers than you have CPU cores. A 4-core machine can handle 5-8 workers without meaningful CPU contention.

**Port range.** Each worker reserves `ports.ports_per_worker` ports (default 10). With the default range of 49152-65535, port exhaustion is not a practical concern below 1000 workers.

```yaml
# Conservative: low-resource machine or rate-limited API
workers:
  max_concurrent: 3

# Aggressive: high-end machine with no rate limits
workers:
  max_concurrent: 10
```

---

## Timeouts

### Per-Task Timeout

`workers.timeout_minutes` sets how long a single task can run before forced termination. The default of 60 minutes is appropriate for most tasks.

| Scenario | Recommended Value | Reasoning |
|----------|------------------:|-----------|
| Simple code generation | 15-30 min | Tasks should not take longer than this |
| Complex refactoring | 45-60 min | Allows for multiple iterations |
| Large file generation | 60-90 min | May involve extensive context reading |
| Integration tasks | 90-120 min | Cross-module coordination is slower |

Tasks that hit the timeout are marked as failed and may be retried (see Retry Configuration below).

### Quality Gate Timeouts

Each gate has its own `timeout` field measured in seconds. Adjust these based on your test suite and tooling speed.

| Gate Type | Recommended Timeout | Notes |
|-----------|--------------------:|-------|
| Lint (ruff, eslint) | 120-300s | Fast for most projects |
| Type check (mypy, tsc) | 300-600s | Scales with codebase size |
| Unit tests | 300-900s | Depends on test count and speed |
| Integration tests | 600-1800s | Include service startup time |
| Security scan | 120-300s | Static analysis is fast |

```yaml
quality_gates:
  - name: test
    command: pytest tests/unit/ -x --timeout=30
    required: true
    timeout: 900  # 15 minutes for a large test suite
```

---

## Retry Configuration

When a task fails (verification command returns non-zero, or worker crashes), ZERG can retry it automatically.

| Option | Effect | Trade-off |
|--------|--------|-----------|
| `retry_attempts: 0` | No retries; fail fast | Fastest feedback, may lose recoverable tasks |
| `retry_attempts: 1` | One retry | Good for transient failures (network, timeout) |
| `retry_attempts: 2` | Two retries (default) | Balanced reliability |
| `retry_attempts: 3-5` | Multiple retries | Use for flaky environments; wastes time on persistent failures |

### Backoff Strategy

The `backoff_strategy` controls the delay between retries.

| Strategy | Behavior | Best For |
|----------|----------|----------|
| `exponential` | 30s, 60s, 120s, ... | Rate limit recovery, transient API errors |
| `linear` | 30s, 60s, 90s, ... | Predictable timing |
| `fixed` | 30s, 30s, 30s, ... | When you know the exact recovery time |

```yaml
# For flaky CI environments
workers:
  retry_attempts: 3
  backoff_strategy: exponential
  backoff_base_seconds: 15
  backoff_max_seconds: 120
```

---

## Context Threshold

`workers.context_threshold_percent` determines when a worker should checkpoint and exit before running out of context window. The context tracker uses heuristics based on files read, tasks executed, tool calls, and elapsed time.

| Value | Behavior |
|------:|----------|
| `50` | Very conservative; workers checkpoint early, more restarts |
| `70` | Default; good balance between context usage and restart overhead |
| `85` | Aggressive; workers use more context, risk truncation |
| `95` | Dangerous; workers may hit hard context limits |

Lower values cause more worker restarts but provide a safety margin. Higher values reduce restarts but risk context overflow, which can cause workers to lose track of their task.

```yaml
# For tasks with many large files
workers:
  context_threshold_percent: 60

# For tasks with few, small files
workers:
  context_threshold_percent: 80
```

---

## Resource Limits

Resource settings primarily affect container mode workers. Subprocess workers are limited by the host OS.

### Container Mode Tuning

```yaml
resources:
  container_memory_limit: "4g"   # Docker --memory flag
  container_cpu_limit: 2.0       # Docker --cpus flag
```

**Memory.** Claude Code workers typically use 300-600 MB. The 4 GB default provides generous headroom for npm installs, test execution, and build steps that happen inside the container.

**CPU.** Workers are I/O-bound. Allocating 1-2 CPU cores per container is sufficient. Allocating more does not meaningfully improve performance.

### When to Increase Resources

- Workers run `npm install` or `pip install` with large dependency trees: increase memory to 6-8 GB.
- Workers compile native extensions (C/C++ deps): increase CPU to 4 and memory to 8 GB.
- Workers run heavy test suites inside the container: increase memory and timeout.

### When to Decrease Resources

- Running many workers on a constrained machine: reduce to 2 GB memory and 1.0 CPU per container.
- Tasks are purely code generation with no build steps: 2 GB memory is sufficient.

---

## Quality Gate Ordering

Gates run sequentially in the order they appear in the configuration. Arrange them from fastest to slowest so that cheap checks fail early.

```yaml
# Recommended ordering: fast checks first
quality_gates:
  - name: lint           # ~5 seconds
    command: ruff check .
    required: true
    timeout: 120

  - name: typecheck      # ~15 seconds
    command: mypy . --strict
    required: true
    timeout: 300

  - name: test           # ~60 seconds
    command: pytest tests/unit/ -x
    required: true
    timeout: 600

  - name: coverage       # ~90 seconds (superset of test)
    command: pytest tests/unit/ --cov=src --cov-fail-under=80
    required: false
    timeout: 600

  - name: security       # ~10 seconds
    command: ruff check . --select S
    required: false
    timeout: 120
```

Setting a gate as `required: false` means its failure produces a warning in `/zerg:status` output but does not block the merge or prevent the next level from starting.

---

## Launcher Type

| Launcher | Pros | Cons | Use When |
|----------|------|------|----------|
| `subprocess` | Fast startup, no Docker needed, shares host tools | No isolation, uses host filesystem | Local development, trusted code |
| `container` | Full isolation, reproducible environment | Slower startup, requires Docker, image setup | CI/CD, untrusted code, team environments |

```yaml
# Local development
workers:
  launcher_type: subprocess

# CI/CD pipeline
workers:
  launcher_type: container
```

Container mode supports two authentication methods. OAuth mounts `~/.claude` into the container. API key mode passes `ANTHROPIC_API_KEY` as an environment variable. See [[Configuration]] for the `ANTHROPIC_API_KEY` environment variable.

---

## Logging Tuning

| Level | Output Volume | Use When |
|-------|---------------|----------|
| `debug` | Very high; includes tool call details | Diagnosing worker failures |
| `info` | Moderate; task progress and gate results | Normal operation |
| `warning` | Low; only problems | Production/CI where log volume matters |
| `error` | Minimal; only failures | Quiet operation |

```yaml
# Debugging a failing task
logging:
  level: debug
  retain_days: 30

# Normal operation
logging:
  level: info
  retain_days: 7
```

---

## Profiles by Use Case

### Solo Developer, Local Machine

```yaml
workers:
  max_concurrent: 3
  timeout_minutes: 45
  retry_attempts: 1
  launcher_type: subprocess
resources:
  cpu_cores: 2
  memory_gb: 4
logging:
  level: info
  retain_days: 3
```

### Team CI/CD Pipeline

```yaml
workers:
  max_concurrent: 8
  timeout_minutes: 60
  retry_attempts: 2
  launcher_type: container
  backoff_strategy: exponential
resources:
  container_memory_limit: "6g"
  container_cpu_limit: 2.0
security:
  level: strict
  container_readonly: true
logging:
  level: warning
  retain_days: 14
```

### Large Monorepo

```yaml
workers:
  max_concurrent: 12
  timeout_minutes: 90
  retry_attempts: 3
  context_threshold_percent: 60
quality_gates:
  - name: lint
    command: ruff check .
    required: true
    timeout: 600
  - name: test
    command: pytest tests/unit/ -x --timeout=60
    required: true
    timeout: 1200
plugins:
  context_engineering:
    task_context_budget_tokens: 6000
```

---

## See Also

- [[Configuration]] -- Full option reference
- [[Plugin System]] -- Adding custom quality gates
- [[Context Engineering]] -- Token budget tuning
