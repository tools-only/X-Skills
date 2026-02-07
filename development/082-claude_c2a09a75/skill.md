# CLAUDE.md - Synkro Development Guide

## ABSOLUTE RULES — DO NOT BREAK THESE
- **NEVER bump the version in pyproject.toml unless the user explicitly tells you to**
- **NEVER publish to PyPI unless the user explicitly tells you to**
- Only bump version for core library changes (synkro/ package code), NOT for notebooks, examples, docs, or config
- If unsure whether to bump, ASK — do not assume

## Quick Reference

```bash
# BEFORE pushing to PyPI or GitHub - ALWAYS run:
pre-commit run --all-files

# Full publish flow:
1. Bump version in pyproject.toml
2. pre-commit run --all-files
3. python -m build
4. python -m twine upload dist/synkro-X.Y.Z*
5. git add -A && git commit -m "..." && git push origin main
```

## Project Overview

Synkro turns unstructured policy documents into LLM training data through a 4-stage pipeline:
1. **Extract** - Parse rules from policy text into structured LogicMap
2. **Generate** - Create golden scenarios (positive, negative, edge cases)
3. **Synthesize** - Generate multi-turn conversation traces
4. **Verify** - Grade traces for quality and policy compliance

## Architecture

```
synkro/
├── pipeline/runner.py      # Main orchestrator - START HERE to understand flow
├── interactive/
│   ├── live_display.py     # Terminal UI (Rich Live) - most UI bugs are here
│   ├── hitl_session.py     # HITL state management
│   └── rich_ui.py          # Static Rich components
├── reporting.py            # Connects pipeline → UI (callbacks)
├── types/
│   ├── logic_map.py        # Rule, LogicMap, GoldenScenario models
│   └── coverage.py         # CoverageReport model
└── core/
    ├── dataset.py          # Dataset model and serialization
    └── checkpoint.py       # Resume interrupted generations
```

## Critical: Live Display System

The UI uses Rich's `Live` component. **Most UI bugs stem from misunderstanding this.**

### How It Works

```python
# In LiveProgressDisplay.start():
self._live = Live(
    self._render,      # Pass CALLABLE, not result
    refresh_per_second=10,
    transient=True,    # Replace in place, don't stack
)
```

### The Three Rules

1. **Never `console.print()` while Live is running** → causes panel stacking
2. **Pass callable to Live, not result** → `self._render` not `self._render()`
3. **Use `refresh()` not `update()`** → `update()` replaces the callable

### Common Bugs and Fixes

| Symptom | Cause | Fix |
|---------|-------|-----|
| Panels stacking vertically | `console.print()` during Live | Update state instead, or return no-op from spinner() |
| Spinner frozen / not animating | Passed `self._render()` to Live | Pass `self._render` (no parens) |
| Time not updating | `_refresh()` using `update()` | Change to `self._live.refresh()` |
| HITL panels stacking | Each render adds new panel | Clear screen before render: `console.clear()` |

### Key Methods

```python
# Start/stop the live display
display.start(model="gemini-2.5-flash")
display.stop()

# Update state (triggers refresh)
display.update_phase("Generating Traces")
display.update_progress(5, 10)
display.add_event("RULES: Extracted 14 rules")

# For HITL - unified clear + render + input
feedback = display.hitl_get_input(logic_map, scenarios, coverage, turns)
```

## Reporting System

`RichReporter` bridges the pipeline and UI via callbacks:

```python
# In pipeline/runner.py:
self.reporter.on_logic_map_complete(logic_map)  # Updates UI
self.reporter.on_response_progress(5, 10)        # Updates progress

# The reporter updates LiveProgressDisplay state internally
```

**Important**: `reporter.spinner()` returns no-op when Live is active to prevent stacking.

## Testing

```bash
# Run all tests (some need API keys)
pytest tests/ -v

# Run without API-dependent tests
pytest tests/test_imports.py tests/test_types.py -v

# Quick smoke test
python -c "from synkro import Generator; print('OK')"
```

## Do's and Don'ts

### Do
- Run `pre-commit run --all-files` before every publish
- Update `DisplayState` fields, not print directly
- Use `hitl_get_input()` for HITL interaction (clears screen)
- Add events via `display.add_event()` for the activity log
- Check `if self._live is not None` before assuming Live is running

### Don't
- Don't `console.print()` while Live display is active
- Don't use `self._live.update(self._render())` - use `refresh()`
- Don't create Rich `Status`/spinners while Live is running
- Don't forget to bump version before publishing
- Don't push without running pre-commit

## Code Patterns

### Adding a new phase to the UI

```python
# In reporting.py callback:
def on_my_new_phase(self, data) -> None:
    self._display.update_phase("My Phase")
    self._display.add_event(f"EVENT: {data}")
    self._display.set_my_data(data)  # If needed
```

### Adding state to DisplayState

```python
# In live_display.py:
@dataclass
class DisplayState:
    # ... existing fields ...
    my_new_field: int = 0

# Add setter method:
def set_my_data(self, data) -> None:
    self._state.my_new_field = data
    self._refresh()
```

### Conditional behavior when Live is active

```python
if self._display._live is not None:
    # Live is running - don't print, update state instead
    self._display.update_phase("Working")
    return _NoOpContextManager()
else:
    # Live not running - safe to print
    return Status("Working...", console=self.console)
```

## Colab Notebooks

**CRITICAL: Shell commands in Colab use `%` magic, NOT `!`**

```python
# CORRECT - use % for shell commands
%pip install -q synkro
%system curl -fsSL https://ollama.ai/install.sh | sh
%system ollama pull model-name
%system head -c 1000 file.txt

# WRONG - do NOT use !
!pip install synkro  # WRONG
!curl ...            # WRONG
!ollama pull ...     # WRONG
```

## Environment

- **Python**: 3.10+ required (uses `|` union syntax, match statements)
- **Build**: Hatchling (`python -m build`)
- **Lint/Format**: Ruff (via pre-commit)
- **Publish**: Twine to PyPI
