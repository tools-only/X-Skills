# Research – Tools Layer Cross-Layer Import Analysis
**Date:** 2026-01-28
**Owner:** agent
**Phase:** Research

## Goal
Full audit of tools layer imports from configuration, exceptions, constants, and types. Map which files, what logic, and why.

---

## Findings by Layer

### 1. Tools → Configuration (5 files)

| Tool File | Imports From | Symbols | Purpose |
|-----------|--------------|---------|---------|
| `tools/ignore.py:11-14` | `configuration.ignore_patterns` | `DEFAULT_IGNORE_PATTERNS`, `DEFAULT_EXCLUDE_DIRS` | Baseline gitignore rules for glob/grep |
| `tools/read_file.py:6` | `configuration.limits` | `get_read_limit()`, `get_max_line_length()` | Line count and width limits |
| `tools/grep.py:23` | `configuration.defaults` | `DEFAULT_USER_CONFIG` | Ripgrep timeout/metrics config |
| `tools/bash.py:11` | `configuration.limits` | `get_command_limit()` | Command output truncation limit |
| `tools/lsp/diagnostics.py:9` | `configuration.user_config` | `load_config()` | LSP enable flag and timeout |

**Anomaly:** `grep.py` imports `DEFAULT_USER_CONFIG` directly instead of using `load_config()`. User config overrides for ripgrep settings are ignored.

---

### 2. Tools → Exceptions (4 files)

| Tool File | Imports From | Exceptions | Usage |
|-----------|--------------|------------|-------|
| `tools/read_file.py:12` | `exceptions` | `ToolExecutionError` | File too large (>100KB) |
| `tools/decorators.py:17` | `exceptions` | `FileOperationError`, `ToolExecutionError` | Permission/decode/OS errors, catch-all wrapper |
| `tools/grep.py:24` | `exceptions` | `TooBroadPatternError`, `ToolExecutionError` | 3-second timeout, unknown search type |
| `tools/parsing/command_parser.py:13` | `exceptions` | `ValidationError` | Invalid JSON, wrong arg type |

**Exception Hierarchy:**
```
TunaCodeError (base)
├── ToolExecutionError
│   └── TooBroadPatternError (grep-specific)
├── FileOperationError
└── ValidationError
```

**Pattern:** Decorators wrap raw exceptions into structured errors with `tool_name`, `message`, `original_error`.

---

### 3. Tools → Constants (3 files)

| Tool File | Imports From | Symbols | Purpose |
|-----------|--------------|---------|---------|
| `tools/read_file.py:7-11` | `constants` | `ERROR_FILE_TOO_LARGE`, `MAX_FILE_SIZE`, `MSG_FILE_SIZE_LIMIT` | File size validation + error messages |
| `tools/bash.py:12-17` | `constants` | `CMD_OUTPUT_TRUNCATED`, `COMMAND_OUTPUT_END_SIZE`, `COMMAND_OUTPUT_START_INDEX`, `COMMAND_OUTPUT_THRESHOLD` | Output truncation strategy |
| `tools/parsing/command_parser.py:8-12` | `constants` | `JSON_PARSE_BASE_DELAY`, `JSON_PARSE_MAX_DELAY`, `JSON_PARSE_MAX_RETRIES` | JSON retry configuration |

**Values:**
| Constant | Value | Used In |
|----------|-------|---------|
| `MAX_FILE_SIZE` | 102400 (100KB) | read_file.py:51 |
| `COMMAND_OUTPUT_THRESHOLD` | 3500 | bash.py:219 |
| `COMMAND_OUTPUT_START_INDEX` | 2500 | bash.py:218,222 |
| `COMMAND_OUTPUT_END_SIZE` | 1000 | bash.py:220 |
| `JSON_PARSE_MAX_RETRIES` | 10 | command_parser.py:38 |
| `JSON_PARSE_BASE_DELAY` | 0.1s | command_parser.py:39 |
| `JSON_PARSE_MAX_DELAY` | 5.0s | command_parser.py:40 |

---

### 4. Tools → Types (1 file)

| Tool File | Imports From | Symbols | Purpose |
|-----------|--------------|---------|---------|
| `tools/parsing/command_parser.py:14` | `types` | `ToolArgs` | Return type annotation (`dict[str, Any]`) |

---

## Complete Dependency Map

```
tools/
├── ignore.py
│   └── configuration.ignore_patterns (DEFAULT_IGNORE_PATTERNS, DEFAULT_EXCLUDE_DIRS)
│
├── read_file.py
│   ├── configuration.limits (get_read_limit, get_max_line_length)
│   ├── constants (ERROR_FILE_TOO_LARGE, MAX_FILE_SIZE, MSG_FILE_SIZE_LIMIT)
│   └── exceptions (ToolExecutionError)
│
├── grep.py
│   ├── configuration.defaults (DEFAULT_USER_CONFIG)
│   └── exceptions (TooBroadPatternError, ToolExecutionError)
│
├── bash.py
│   ├── configuration.limits (get_command_limit)
│   └── constants (CMD_OUTPUT_TRUNCATED, COMMAND_OUTPUT_*)
│
├── decorators.py
│   └── exceptions (FileOperationError, ToolExecutionError)
│
├── lsp/diagnostics.py
│   └── configuration.user_config (load_config)
│
└── parsing/command_parser.py
    ├── constants (JSON_PARSE_*)
    ├── exceptions (ValidationError)
    └── types (ToolArgs)
```

---

## Dependency Chain Analysis

### Chain 1: tools → configuration → configuration
```
tools/grep.py
  → configuration.defaults (DEFAULT_USER_CONFIG)
    → (self-contained, no further deps)

tools/bash.py
tools/read_file.py
  → configuration.limits
    → configuration.user_config (load_config)
    → constants (DEFAULT_READ_LIMIT, etc.)
```

### Chain 2: tools → exceptions (flat)
```
tools/*
  → exceptions (TunaCodeError hierarchy)
    → (no further deps, exceptions.py is root-level)
```

### Chain 3: tools → constants (flat)
```
tools/*
  → constants
    → (no further deps, constants.py is root-level)
```

### Chain 4: tools → types (flat)
```
tools/parsing/command_parser.py
  → types
    → (types imports from typing only)
```

---

## Potential Issues

### Issue 1: grep.py Bypasses User Config
**File:** `tools/grep.py:49-51`
```python
def _load_ripgrep_config(self) -> dict[str, Any]:
    return DEFAULT_USER_CONFIG["settings"]["ripgrep"]
```
User settings in `~/.tunacode/tunacode.json` are ignored for ripgrep timeout/metrics.

### Issue 2: limits.py Has Delayed Import
**File:** `configuration/limits.py:23-24`
```python
def _load_settings() -> dict:
    from tunacode.configuration.user_config import load_config  # late import
```
Indicates circular import concern between limits ↔ user_config.

### Issue 3: No Cache Invalidation
**File:** `configuration/limits.py:32-34`
```python
def clear_cache() -> None:
    _load_settings.cache_clear()
```
`clear_cache()` exists but is never called. Config changes during session use stale cached values.

---

## Summary Statistics

| Layer | Files Importing | Total Symbols | Notes |
|-------|-----------------|---------------|-------|
| configuration | 5 | 9 | Most diverse: patterns, limits, defaults, user_config |
| exceptions | 4 | 4 | Clean hierarchy, used for error wrapping |
| constants | 3 | 10 | Operational limits and error templates |
| types | 1 | 1 | Only ToolArgs type alias |

**Total cross-layer imports in tools/:** 13 files touching 24 symbols from 4 layers.

---

## References

- `src/tunacode/tools/ignore.py`
- `src/tunacode/tools/read_file.py`
- `src/tunacode/tools/grep.py`
- `src/tunacode/tools/bash.py`
- `src/tunacode/tools/decorators.py`
- `src/tunacode/tools/lsp/diagnostics.py`
- `src/tunacode/tools/parsing/command_parser.py`
- `src/tunacode/configuration/ignore_patterns.py`
- `src/tunacode/configuration/limits.py`
- `src/tunacode/configuration/defaults.py`
- `src/tunacode/configuration/user_config.py`
- `src/tunacode/exceptions.py`
- `src/tunacode/constants.py`
- `src/tunacode/types/__init__.py`
