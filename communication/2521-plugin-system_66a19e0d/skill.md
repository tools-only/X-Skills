# Plugin System

ZERG supports plugins for extending quality gates, reacting to lifecycle events, and adding custom worker launcher backends. Plugins are additive only -- they observe and extend, but cannot mutate orchestrator state. Plugin failures are isolated and never crash the orchestrator.

For the full API reference, see [[Plugin API Reference]]. For configuration options, see [[Configuration]].

---

## Plugin Types

ZERG provides four plugin extension points:

| Plugin Type | Purpose | Blocks Execution | Example Use Case |
|-------------|---------|:-----------------:|------------------|
| Quality Gate | Validate code after merges | Yes (if `required`) | Security scan, complexity check, license audit |
| Lifecycle Hook | React to ZERG events | No | Slack notifications, metrics reporting, audit logs |
| Launcher | Custom worker environments | No | Kubernetes pods, SSH clusters, cloud VMs |
| Context | Custom context injection | No | Project-specific context, custom summarization |

---

## Quality Gate Plugins

Quality gates validate code after each level merge. They run sequentially after the built-in gates defined in `quality_gates` of the config.

### YAML Configuration (Simple)

For gates that are shell commands, define them directly in `.zerg/config.yaml`:

```yaml
plugins:
  quality_gates:
    - name: security-scan
      command: bandit -r src/ --severity medium
      required: false
      timeout: 300

    - name: complexity-check
      command: radon cc src/ --min B
      required: true
      timeout: 120

    - name: license-audit
      command: pip-licenses --fail-on "GPL-3.0"
      required: false
      timeout: 60
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `name` | string | required | Unique gate identifier |
| `command` | string | required | Shell command to execute |
| `required` | bool | `false` | If true, failure blocks the merge |
| `timeout` | int | `300` | Max seconds before timeout (1-3600) |

### Python Plugin (Advanced)

For gates that need custom logic, API calls, or multi-step validation, implement the `QualityGatePlugin` abstract base class. See [[Plugin API Reference]] for the full class definition.

```python
from zerg.plugins import QualityGatePlugin, GateContext
from zerg.types import GateRunResult
from zerg.constants import GateResult

class PerformanceBenchmark(QualityGatePlugin):
    @property
    def name(self) -> str:
        return "perf-benchmark"

    def run(self, ctx: GateContext) -> GateRunResult:
        # Run benchmark and compare against baseline
        # ...
        return GateRunResult(
            gate_name=self.name,
            result=GateResult.PASS,
            command="bench --compare baseline",
            exit_code=0,
            stdout="All benchmarks within threshold",
            stderr="",
        )
```

### Execution Order

```
Built-in gates (from quality_gates config)
  lint -> typecheck -> test -> coverage -> security
        |
        v
Plugin gates (from plugins.quality_gates config)
  security-scan -> complexity-check -> license-audit
        |
        v
Python plugin gates (from entry points)
  perf-benchmark -> ...
        |
        v
Merge completes (if all required gates passed)
```

---

## Lifecycle Hook Plugins

Hooks observe events during the ZERG lifecycle without blocking execution. They are useful for notifications, logging, and metrics collection.

### Available Events

| Event Type | Emitted When | Data Payload |
|------------|-------------|--------------|
| `task_started` | Worker begins task execution | `task_id`, `worker_id`, `level` |
| `task_completed` | Task verification passes | `task_id`, `worker_id`, `duration`, `output` |
| `level_complete` | All tasks in a level finish | `level`, `task_count`, `elapsed_time` |
| `merge_complete` | Level branches merged | `level`, `branch_count`, `conflicts` |
| `worker_spawned` | New worker process starts | `worker_id`, `port`, `mode` |
| `quality_gate_run` | Quality gate executes | `gate_name`, `result`, `level` |
| `rush_started` | `/zerg:rush` begins | `feature`, `workers`, `config` |
| `rush_finished` | All levels complete | `feature`, `total_time`, `tasks_completed` |

### YAML Configuration (Simple)

For hooks that execute shell commands:

```yaml
plugins:
  hooks:
    - event: task_completed
      command: echo "Task {task_id} completed at $(date)"
      timeout: 60

    - event: level_complete
      command: ./scripts/notify-slack.sh "Level {level} done"
      timeout: 120

    - event: merge_complete
      command: |
        python scripts/generate_report.py \
          --level {level} \
          --feature {feature}
      timeout: 180

    - event: rush_finished
      command: ./scripts/send-summary.sh "{feature}" "{total_time}"
      timeout: 60
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `event` | string | required | Event type from the table above |
| `command` | string | required | Shell command (parsed with `shlex.split`, no shell=True) |
| `timeout` | int | `60` | Max seconds before timeout (1-600) |

**Variable substitution.** Curly-brace placeholders like `{level}`, `{feature}`, `{task_id}`, and `{worker_id}` are replaced with values from the event data payload.

### Python Plugin (Advanced)

For hooks with custom logic, implement the `LifecycleHookPlugin` abstract base class:

```python
from zerg.plugins import LifecycleHookPlugin, LifecycleEvent

class SlackNotifier(LifecycleHookPlugin):
    @property
    def name(self) -> str:
        return "slack-notifier"

    def on_event(self, event: LifecycleEvent) -> None:
        if event.event_type == "rush_finished":
            send_slack_message(
                f"Feature {event.data['feature']} built in "
                f"{event.data['total_time']}s"
            )
```

### Safety Guarantees

- Hooks never block orchestrator execution.
- Each hook invocation is wrapped in its own try/except.
- One failing hook does not prevent other hooks from running.
- Hooks have no access to mutable orchestrator state.

---

## Launcher Plugins

Launcher plugins add custom worker execution environments beyond the built-in `subprocess` and `container` modes.

### YAML Configuration

```yaml
plugins:
  launchers:
    - name: kubernetes
      entry_point: my_pkg.launchers:K8sLauncher
```

| Field | Type | Description |
|-------|------|-------------|
| `name` | string | Launcher name referenced in `workers.launcher_type` |
| `entry_point` | string | Python entry point in `module:Class` format |

After registering a launcher plugin, use it by setting:

```yaml
workers:
  launcher_type: kubernetes
```

### Python Plugin (Advanced)

```python
from zerg.plugins import LauncherPlugin

class K8sLauncherPlugin(LauncherPlugin):
    @property
    def name(self) -> str:
        return "kubernetes"

    def create_launcher(self, config):
        return K8sWorkerLauncher(config)
```

The returned object must implement the `WorkerLauncher` abstract base class with `launch()`, `wait()`, and `cleanup()` methods. See [[Plugin API Reference]] for details.

---

## Context Plugins

Context plugins customize the context injected into worker prompts. The built-in `ContextEngineeringPlugin` handles security rule filtering and spec excerpts. Custom context plugins can add project-specific context.

See [[Context Engineering]] for the built-in plugin and [[Plugin API Reference]] for the abstract base class.

---

## Plugin Registration

Plugins are discovered and registered through three mechanisms, in order of priority:

### 1. YAML Configuration

The simplest approach. Define hooks and gates directly in `.zerg/config.yaml`. Suitable for shell commands and simple automation.

### 2. Python Entry Points

For installable plugins distributed as Python packages. Register plugins via `pyproject.toml`:

```toml
[project.entry-points."zerg.plugins"]
my-gate = "my_package.gates:MyGatePlugin"
my-hook = "my_package.hooks:MyHookPlugin"
my-launcher = "my_package.launchers:MyLauncherPlugin"
```

ZERG discovers entry points in the `zerg.plugins` group at startup and registers them based on which abstract base class they implement.

### 3. Direct Registration

For programmatic use, register plugins directly with the `PluginRegistry`:

```python
from zerg.plugins import PluginRegistry

registry = PluginRegistry()
registry.register_gate(MyGatePlugin(), required=True)
registry.register_hook("task_completed", my_callback)
registry.register_launcher(MyLauncherPlugin())
```

---

## Task System Integration

All plugin operations integrate with the Claude Code Task system per ZERG conventions.

| Operation | Task Subject |
|-----------|-------------|
| Plugin registration | `[Plugins] Register {name}` |
| Quality gate execution | `[Gate] Run {gate_name}` |
| Lifecycle hook invocation | `[Hook] Process {event_type}` |
| Custom launcher execution | `[Launcher] Spawn {worker_id}` |

Plugin state is tracked in both the Task system (authoritative) and `.zerg/state/plugins.json` (supplementary). If they disagree, the Task system wins.

---

## See Also

- [[Plugin API Reference]] -- Full class definitions and method signatures
- [[Configuration]] -- All configuration options
- [[Context Engineering]] -- Built-in context engineering plugin
