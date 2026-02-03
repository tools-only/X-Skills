# Research – Utils Cross-Layer Import Analysis
**Date:** 2026-01-28
**Owner:** agent
**Phase:** Research

## Goal
Map out the three cross-layer imports involving utils, configuration, and types modules. Document the actual import direction, files involved, and logic shared—without judgment on correctness.

## Findings

### Import 1: utils.system.gitignore ← configuration.ignore_patterns

**Direction:** utils imports FROM configuration (not the reverse)

**File:** `src/tunacode/utils/system/gitignore.py:10-14`

```python
from tunacode.configuration.ignore_patterns import (
    DEFAULT_IGNORE_PATTERNS,
    GIT_DIR_PATTERN,
    is_ignored,
)
```

**What's imported:**
| Symbol | Source Location | Type |
|--------|-----------------|------|
| `DEFAULT_IGNORE_PATTERNS` | `ignore_patterns.py:15-68` | Tuple of 53 ignore patterns |
| `GIT_DIR_PATTERN` | `ignore_patterns.py:12` | String constant `".git/"` |
| `is_ignored()` | `ignore_patterns.py:84-126` | Function (fnmatch logic) |

**Logic/Usage:**
- `DEFAULT_IGNORE_PATTERNS`: Fallback when no `.gitignore` exists (`gitignore.py:51`)
- `GIT_DIR_PATTERN`: Always added to patterns, even if not in `.gitignore` (`gitignore.py:31`)
- `is_ignored()`: Core filtering in `list_cwd()` to check dirs/files against patterns (`gitignore.py:72, 78`)

---

### Import 2: utils.system ← configuration.paths

**Direction:** utils re-exports FROM configuration (not the reverse)

**File:** `src/tunacode/utils/system/__init__.py:3-9`

```python
from tunacode.configuration.paths import (
    check_for_updates,
    cleanup_session,
    get_cwd,
    get_session_dir,
    get_tunacode_home,
)
```

**What's imported:**
| Symbol | Source Location | Type |
|--------|-----------------|------|
| `check_for_updates()` | `paths.py:131-160` | Function (PyPI version check) |
| `cleanup_session()` | `paths.py:81-106` | Function (delete temp files) |
| `get_cwd()` | `paths.py:41-43` | Function (returns cwd) |
| `get_session_dir()` | `paths.py:26-38` | Function (returns session path) |
| `get_tunacode_home()` | `paths.py:13-23` | Function (returns `~/.tunacode`) |

**Logic/Usage:**
- These are pure filesystem/system utilities
- `utils.system` re-exports them via `__all__` (lines 15-23)
- Consumers can import from `tunacode.utils.system` instead of `tunacode.configuration.paths`

**Note:** `configuration/paths.py` itself imports from `configuration.settings` (line 9) for version check.

---

### Import 3: utils.messaging.adapter ← types.canonical

**Direction:** utils imports FROM types

**File:** `src/tunacode/utils/messaging/adapter.py:18-28`

```python
from tunacode.types.canonical import (
    CanonicalMessage,
    CanonicalPart,
    MessageRole,
    RetryPromptPart,
    SystemPromptPart,
    TextPart,
    ThoughtPart,
    ToolCallPart,
    ToolReturnPart,
)
```

**What's imported:**
| Symbol | Type |
|--------|------|
| `CanonicalMessage` | Frozen dataclass (message container) |
| `CanonicalPart` | Union type alias |
| `MessageRole` | Enum (USER/ASSISTANT/TOOL/SYSTEM) |
| `RetryPromptPart` | Frozen dataclass |
| `SystemPromptPart` | Frozen dataclass |
| `TextPart` | Frozen dataclass |
| `ThoughtPart` | Frozen dataclass |
| `ToolCallPart` | Frozen dataclass |
| `ToolReturnPart` | Frozen dataclass |

**Logic/Usage:**
- `_determine_role()` returns `MessageRole` enum values (`adapter.py:119-138`)
- `_convert_part_to_canonical()` instantiates part dataclasses (`adapter.py:75-116`)
- `to_canonical()` creates `CanonicalMessage` instances (`adapter.py:146-203`)
- `from_canonical()` pattern matches on part types via `isinstance` (`adapter.py:215-268`)
- `get_content()` checks `isinstance(message, CanonicalMessage)` (`adapter.py:280-288`)

---

## Summary Table

| Import | Direction | Importing Module | Imported From |
|--------|-----------|------------------|---------------|
| 1 | utils ← config | `utils/system/gitignore.py` | `configuration/ignore_patterns.py` |
| 2 | utils ← config | `utils/system/__init__.py` | `configuration/paths.py` |
| 3 | utils ← types | `utils/messaging/adapter.py` | `types/canonical.py` |

## Dependency Chain

```
utils.system.gitignore → configuration.ignore_patterns
utils.system            → configuration.paths → configuration.settings
utils.messaging.adapter → types.canonical
```

## References
- `src/tunacode/utils/system/gitignore.py`
- `src/tunacode/utils/system/__init__.py`
- `src/tunacode/utils/messaging/adapter.py`
- `src/tunacode/configuration/ignore_patterns.py`
- `src/tunacode/configuration/paths.py`
- `src/tunacode/types/canonical.py`
