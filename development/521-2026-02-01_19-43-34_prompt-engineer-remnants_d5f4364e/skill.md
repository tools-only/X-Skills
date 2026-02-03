# Research - Prompt Engineer Remnants Removal

**Date:** 2026-02-01
**Owner:** agent
**Phase:** Research
**Git Commit:** 35aa99e2

## Goal

Identify all remnants of the "prompt engineer" feature that can be removed now that the codebase has returned to a simple prompt system. Exclude toolprompts from consideration.

## Findings

### Dead Code (Safe to Remove)

| File | Location | Description |
|------|----------|-------------|
| `src/tunacode/core/prompting/prompting_engine.py` | Lines 82-98 | `compose_prompt()` function - never called, not exported, not used |

### Over-Engineered Code (Consider Simplifying)

| File | Location | Description |
|------|----------|-------------|
| `src/tunacode/core/prompting/prompting_engine.py` | Lines 33-40 | `register()` method - tested but never used in production |
| `src/tunacode/core/prompting/prompting_engine.py` | Lines 65-74 | Singleton pattern (`_engine`, `get_prompting_engine()`) - only used internally |
| `src/tunacode/core/prompting/prompting_engine.py` | Lines 10-62 | `PromptingEngine` class - could be replaced with a simple function |

### What IS Used (Do Not Remove)

- `resolve_prompt()` - Called from `agent_config.py:225` to resolve `{{CWD}}`, `{{OS}}`, `{{DATE}}` placeholders
- Built-in placeholder providers (CWD, OS, DATE) - Used in `system_prompt.md:124-127`

### Test File

| File | Notes |
|------|-------|
| `tests/unit/core/test_prompting_engine.py` | Tests for unused features (custom provider, singleton) will need removal/update |

### Documentation

| File | Notes |
|------|-------|
| `docs/codebase-map/modules/core-prompting.md` | Describes unused extension points (lines 56-66) |
| `docs/codebase-map/modules/00-overview.md` | Line 18 mentions "prompt engineering" |

## Recommended Removal Strategy

### Option A: Minimal (Remove Dead Code Only)

1. Delete `compose_prompt()` function (lines 82-98)
2. Keep everything else as-is

### Option B: Simplify (Flatten the Abstraction)

Replace the entire `PromptingEngine` class with a simple function:

```python
"""Prompting engine for resolving {{placeholder}} syntax in prompts."""

import os
import platform
import re
from datetime import datetime

_PLACEHOLDER_PATTERN = re.compile(r"\{\{(.+?)\}\}")

_PROVIDERS = {
    "CWD": os.getcwd,
    "OS": platform.system,
    "DATE": lambda: datetime.now().isoformat(),
}

def resolve_prompt(template: str) -> str:
    """Resolve all placeholders in the template."""
    if not template:
        return template

    def replace_match(match: re.Match[str]) -> str:
        name = match.group(1).strip()
        provider = _PROVIDERS.get(name)
        return provider() if provider else match.group(0)

    return _PLACEHOLDER_PATTERN.sub(replace_match, template)
```

This removes:
- `PromptingEngine` class
- `compose_prompt()` function
- `register()` method
- Singleton pattern
- ~50 lines of code

### Option C: Complete Removal (Inline Everything)

If placeholders are rarely changed, inline the resolution directly in `agent_config.py:load_system_prompt()`:

```python
def load_system_prompt(base_path: Path, model: str | None = None) -> str:
    prompt_file = base_path / "prompts" / "system_prompt.md"
    if not prompt_file.exists():
        raise FileNotFoundError(f"Required prompt file not found: {prompt_file}")

    content = prompt_file.read_text(encoding="utf-8")

    # Inline placeholder resolution
    replacements = {
        "{{CWD}}": os.getcwd(),
        "{{OS}}": platform.system(),
        "{{DATE}}": datetime.now().isoformat(),
    }
    for placeholder, value in replacements.items():
        content = content.replace(placeholder, value)

    return content
```

This would allow deleting the entire `src/tunacode/core/prompting/` directory.

## Files Affected by Each Option

### Option A
- `src/tunacode/core/prompting/prompting_engine.py` - Remove lines 82-98

### Option B
- `src/tunacode/core/prompting/prompting_engine.py` - Replace entirely
- `tests/unit/core/test_prompting_engine.py` - Update tests (remove class-based tests)

### Option C
- `src/tunacode/core/prompting/` - Delete entire directory
- `src/tunacode/core/agents/agent_components/agent_config.py` - Update import, inline logic
- `tests/unit/core/test_prompting_engine.py` - Delete file or move tests
- `docs/codebase-map/modules/core-prompting.md` - Delete file

## Knowledge Gaps

- Confirm whether custom placeholder registration was ever planned for future use
- Verify no other files import from `prompting_engine` besides `agent_config.py`

## References

- [`src/tunacode/core/prompting/prompting_engine.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/35aa99e2/src/tunacode/core/prompting/prompting_engine.py)
- [`src/tunacode/core/prompting/__init__.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/35aa99e2/src/tunacode/core/prompting/__init__.py)
- [`tests/unit/core/test_prompting_engine.py`](https://github.com/alchemiststudiosDOTai/tunacode/blob/35aa99e2/tests/unit/core/test_prompting_engine.py)
- [`docs/codebase-map/modules/core-prompting.md`](https://github.com/alchemiststudiosDOTai/tunacode/blob/35aa99e2/docs/codebase-map/modules/core-prompting.md)
