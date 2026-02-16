# Technical Design: launcher-mode-verification

## Metadata
- **Feature**: launcher-mode-verification
- **Status**: APPROVED
- **Created**: 2026-01-31

## 1. Overview

Enhance rush output to show human-friendly launcher mode and worker count at start and in completion summary.

## 2. Changes

### 2.1 `zerg/commands/rush.py` (~lines 159-160, 188-196)

**Start output** — replace:
```python
launcher_name = type(orchestrator.launcher).__name__
console.print(f"Launcher: [bold]{launcher_name}[/bold]")
```
With:
```python
launcher_name = type(orchestrator.launcher).__name__
mode_label = "container (Docker)" if "Container" in launcher_name else "subprocess"
console.print(f"Launcher mode: [bold]{mode_label}[/bold]")
console.print(f"Workers: [bold]{workers}[/bold]")
```

**Completion output** — add mode to summary:
```python
if status["is_complete"]:
    console.print(f"\n[bold green]✓ All tasks complete![/bold green] (mode: {mode_label})")
```

### 2.2 Tests — `tests/unit/test_rush_cmd.py`

Update any assertions that check for the old "Launcher:" output format.

## 3. File Ownership

| File | Task ID | Operation |
|------|---------|-----------|
| zerg/commands/rush.py | TASK-001 | modify |
| tests/unit/test_rush_cmd.py | TASK-002 | modify |

## 4. Verification

```bash
pytest tests/unit/test_rush_cmd.py -v
pytest tests/ -x -q
```
