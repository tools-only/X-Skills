---
title: Utils Limits
path: src/tunacode/configuration/limits.py
type: file
depth: 1
description: Centralized tool output limits with cascading defaults
exports: [get_command_limit, get_max_files_in_dir, get_max_tokens, clear_cache]
seams: [M]
---

# Utils Limits

## Purpose

Centralizes tool output limit configuration with a two-tier precedence system.

## Precedence System

```
explicit setting > standard default
```

Implemented in `_get_limit()`:

```python
def _get_limit(key: str, default: int) -> int:
    settings = _load_settings()

    # Explicit setting wins
    if key in settings:
        return settings[key]

    return default
```

## Key Functions

### get_command_limit()

Max chars from `bash` output.
- Standard: 5000 chars

### get_max_files_in_dir()

Max entries from `list_dir`.
- Standard: 50 files

### get_max_tokens()

Cap on model response length.
- Explicit `max_tokens` setting takes precedence
- Standard mode: returns `None` (unlimited)

### clear_cache()

Clears the `@lru_cache` on settings. Call when config changes at runtime.

## Integration Points

| Consumer | File | Usage |
|----------|------|-------|
| Bash tool | `tools/bash.py` | `get_command_limit()` |
| List dir tool | `tools/list_dir.py` | `get_max_files_in_dir()` |
| Agent config | `core/agents/agent_components/agent_config.py` | `get_max_tokens()` |

## Design Notes

- **Caching**: Uses `@lru_cache(maxsize=1)` for settings load
- **Circular import avoidance**: Imports config loader inside function

## Constants Location

Default values defined in `constants.py`:

```python
MAX_COMMAND_OUTPUT = 5000
MAX_FILES_IN_DIR = 50
```

## Seams (M)

**Modification Points:**
- Add new limit types (follow `_get_limit()` pattern)
- Change default values in `constants.py`
