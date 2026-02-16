# WM-006: Add metrics display to status command

## Overview
Update the status command to display worker and level metrics.

## Files to Modify
- `zerg/commands/status.py`

## Requirements

### New Import
Add with other imports:
```python
from zerg.metrics import MetricsCollector
```

### Helper Function for Duration Formatting

Add this helper function:

```python
def format_duration(ms: int | None) -> str:
    """Format milliseconds as human-readable duration.

    Args:
        ms: Duration in milliseconds

    Returns:
        Formatted string like "1.2s", "4m30s", "1h15m"
    """
    if ms is None:
        return "-"

    seconds = ms / 1000

    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        minutes = int(seconds // 60)
        secs = int(seconds % 60)
        return f"{minutes}m{secs}s"
    else:
        hours = int(seconds // 3600)
        minutes = int((seconds % 3600) // 60)
        return f"{hours}h{minutes}m"
```

### New Function: show_worker_metrics

```python
def show_worker_metrics(state: StateManager) -> None:
    """Show worker metrics table.

    Args:
        state: State manager
    """
    console.print("[bold]Worker Metrics:[/bold]")

    # Compute metrics
    collector = MetricsCollector(state)

    table = Table(show_header=True)
    table.add_column("Worker", justify="center")
    table.add_column("Init", justify="right")
    table.add_column("Tasks", justify="center")
    table.add_column("Avg", justify="right")
    table.add_column("Uptime", justify="right")
    table.add_column("Status")

    workers = state.get_all_workers()

    if not workers:
        table.add_row("-", "-", "-", "-", "-", "[dim]No workers[/dim]")
    else:
        for worker_id in sorted(workers.keys()):
            worker = workers[worker_id]
            metrics = collector.compute_worker_metrics(worker_id)

            # Format status
            if worker.status == WorkerStatus.RUNNING:
                status_display = "[green]RUNNING[/green]"
            elif worker.status == WorkerStatus.IDLE:
                status_display = "[dim]IDLE[/dim]"
            elif worker.status == WorkerStatus.CRASHED:
                status_display = "[red]CRASHED[/red]"
            else:
                status_display = worker.status.value.upper()

            table.add_row(
                f"worker-{worker_id}",
                format_duration(metrics.initialization_ms),
                str(metrics.tasks_completed),
                format_duration(int(metrics.avg_task_duration_ms)) if metrics.avg_task_duration_ms else "-",
                format_duration(metrics.uptime_ms),
                status_display,
            )

    console.print(table)
    console.print()
```

### New Function: show_level_metrics

```python
def show_level_metrics(state: StateManager, level_filter: int | None) -> None:
    """Show level metrics with timing.

    Args:
        state: State manager
        level_filter: Level to filter to
    """
    console.print("[bold]Level Metrics:[/bold]")

    collector = MetricsCollector(state)

    table = Table(show_header=True)
    table.add_column("Level", justify="center")
    table.add_column("Name")
    table.add_column("Duration", justify="right")
    table.add_column("Tasks", justify="center")
    table.add_column("Avg", justify="right")
    table.add_column("p50", justify="right")
    table.add_column("p95", justify="right")
    table.add_column("Status")

    levels = state._state.get("levels", {})
    current_level = state.get_current_level()

    # If no level data, show placeholders
    if not levels:
        for i in range(1, 6):
            if level_filter and i != level_filter:
                continue
            status = "[yellow]RUNNING[/yellow]" if i == current_level else "PENDING"
            table.add_row(str(i), f"Level {i}", "-", "-", "-", "-", "-", status)
    else:
        for level_str in sorted(levels.keys(), key=int):
            level_num = int(level_str)
            if level_filter and level_num != level_filter:
                continue

            level_data = levels[level_str]
            metrics = collector.compute_level_metrics(level_num)

            name = level_data.get("name", f"Level {level_num}")
            status_str = level_data.get("status", "pending")

            if status_str == "complete":
                status_display = "[green]✓ DONE[/green]"
            elif status_str == "running":
                status_display = "[yellow]RUNNING[/yellow]"
            else:
                status_display = "PENDING"

            table.add_row(
                level_str,
                name,
                format_duration(metrics.duration_ms),
                f"{metrics.completed_count}/{metrics.task_count}",
                format_duration(int(metrics.avg_task_duration_ms)) if metrics.avg_task_duration_ms else "-",
                format_duration(metrics.p50_duration_ms) if metrics.p50_duration_ms else "-",
                format_duration(metrics.p95_duration_ms) if metrics.p95_duration_ms else "-",
                status_display,
            )

    console.print(table)
    console.print()
```

### Update show_status Function

Replace the calls to `show_level_status` and `show_worker_status` with the new metrics functions:

```python
def show_status(state: StateManager, feature: str, level_filter: int | None) -> None:
    """Show current status.

    Args:
        state: State manager
        feature: Feature name
        level_filter: Level to filter to
    """
    # Header
    console.print()
    console.print(Panel(f"[bold cyan]ZERG Status: {feature}[/bold cyan]"))
    console.print()

    # Get task stats
    all_tasks = state._state.get("tasks", {})
    completed = len(state.get_tasks_by_status(TaskStatus.COMPLETE))
    failed = len(state.get_tasks_by_status(TaskStatus.FAILED))
    total = len(all_tasks) if all_tasks else 0

    # Progress bar
    progress_pct = (completed / total * 100) if total > 0 else 0
    progress_bar = create_progress_bar(progress_pct)
    console.print(f"Progress: {progress_bar} {progress_pct:.0f}% ({completed}/{total} tasks)")
    console.print()

    # Level metrics (replaces show_level_status)
    show_level_metrics(state, level_filter)

    # Worker metrics (replaces show_worker_status)
    show_worker_metrics(state)

    # Compute feature-level metrics for summary
    collector = MetricsCollector(state)
    metrics = collector.compute_feature_metrics()

    # Summary line
    console.print(
        f"Total: {format_duration(metrics.total_duration_ms)} | "
        f"Workers: {metrics.workers_used} | "
        f"Tasks: {metrics.tasks_completed}/{metrics.tasks_total}"
    )
    console.print()

    # Recent events
    show_recent_events(state, limit=5)

    # Error state
    error = state.get_error()
    if error:
        console.print(f"\n[red]Error:[/red] {error}")
```

### Update show_json_status

Include metrics in JSON output:

```python
def show_json_status(state: StateManager, level_filter: int | None) -> None:
    """Output status as JSON.

    Args:
        state: State manager
        level_filter: Level to filter to
    """
    # Compute metrics
    collector = MetricsCollector(state)
    metrics = collector.compute_feature_metrics()

    output = {
        "feature": state.feature,
        "current_level": state.get_current_level(),
        "paused": state.is_paused(),
        "error": state.get_error(),
        "tasks": state._state.get("tasks", {}),
        "workers": {
            str(wid): w.to_dict() for wid, w in state.get_all_workers().items()
        },
        "levels": state._state.get("levels", {}),
        "events": state.get_events(limit=10),
        "metrics": metrics.to_dict(),  # Add metrics
    }

    console.print(json.dumps(output, indent=2, default=str))
```

## Verification

```bash
python -c "from zerg.commands.status import show_status; print('Status import OK')"
```

## Acceptance Criteria
- [ ] format_duration helper function exists
- [ ] show_worker_metrics displays worker timing table
- [ ] show_level_metrics displays level timing with p50/p95
- [ ] show_status calls the new metrics functions
- [ ] JSON output includes metrics section
- [ ] Graceful handling when no metrics available (shows "-")
- [ ] Durations formatted as human-readable (1.2s, 4m30s)
