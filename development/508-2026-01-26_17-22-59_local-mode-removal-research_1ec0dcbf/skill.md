# Research – local_mode Feature Removal

**Date:** 2026-01-26
**Owner:** tunacode
**Phase:** Research
**git_commit:** b6eb3990

## Goal

Document the `local_mode` feature for complete removal. The feature was introduced when the codebase was messy to help users with small context window models. Users can now pass their own `OPENAI_BASE_URL` to use local models without this feature.

## Original Commit

**Commit:** `cb5a4568a812696e7d41796736a498cf77b6ca96`
**Date:** 2026-01-07
**Author:** tunahorse1
**Message:** feat: add local_mode for small context window models

**Original changes:**
- Added `LOCAL_TEMPLATE` with minimal sections (AGENT_ROLE, TOOL_USE, USER_INSTRUCTIONS)
- Reduced tools from 11 to 6 in local mode (no grep, web_fetch, research, todo)
- Used 1-word tool descriptions to save tokens
- Added aggressive pruning thresholds (2k protect vs 40k, 500 min vs 20k)
- Added local limits: read 200 lines, bash 1500 chars, line 500 chars
- Added `local_max_tokens` setting to cap response length
- Made `guide_file` setting work
- Created `CLAUDE_LOCAL.md` with condensed instructions

**Files modified in original commit:**
- `.claude/JOURNAL.md`
- `CLAUDE_LOCAL.md` (deleted - no longer exists)
- `src/tunacode/constants.py`
- `src/tunacode/core/agents/agent_components/agent_config.py`
- `src/tunacode/core/compaction.py`
- `src/tunacode/core/prompting/__init__.py`
- `src/tunacode/core/prompting/templates.py`
- `src/tunacode/tools/bash.py`
- `src/tunacode/tools/read_file.py`
- `tests/test_compaction.py`

## Feature Scope

The `local_mode` feature was a system-wide optimization toggle that affected 6 layers:

| Layer | Standard Mode | Local Mode |
|-------|--------------|------------|
| 1. System prompt | MAIN_TEMPLATE (5 sections) | LOCAL_TEMPLATE (3 sections) |
| 2. Guide file | AGENTS.md (~2k tokens) | CLAUDE_LOCAL.md (~500 tokens) |
| 3. Tool schemas | 11 tools (~1.8k tokens) | 6 tools (~575 tokens) |
| 4. Output limits | read: 2000, bash: 5000 | read: 200, bash: 1500 |
| 5. Response cap | Unlimited | 1000 tokens |
| 6. Pruning | Protect 40k tokens | Protect 2k tokens |

## Findings

### Production Code Files (src/tunacode/)

#### 1. `src/tunacode/constants.py:34-38`

**Lines:**
```python
# Local mode limits (for small context windows)
LOCAL_MAX_COMMAND_OUTPUT = 1500
LOCAL_DEFAULT_READ_LIMIT = 200
LOCAL_MAX_LINE_LENGTH = 500
LOCAL_MAX_FILES_IN_DIR = 20
```

**Action:** DELETE lines 34-38

---

#### 2. `src/tunacode/utils/limits.py` (ENTIRE FILE)

This is the central control point for `local_mode`. The entire file's logic is built around the precedence system: `explicit setting > local_mode default > standard default`.

**Lines to modify:**
- **Lines 3-8:** Delete comment about local_mode precedence
- **Lines 14-15:** Delete `LOCAL_DEFAULT_READ_LIMIT`, `LOCAL_MAX_*` imports
- **Lines 42-56:** Simplify `_get_limit()` - remove local_mode logic
- **Lines 59-61:** DELETE `is_local_mode()` function entirely
- **Lines 64-66:** Update `get_read_limit()` to use simplified logic
- **Lines 69-71:** Update `get_max_line_length()` to use simplified logic
- **Lines 74-76:** Update `get_command_limit()` to use simplified logic
- **Lines 79-81:** Update `get_max_files_in_dir()` to use simplified logic
- **Lines 84-96:** DELETE or simplify `get_max_tokens()` - remove local_mode logic

**Action:** MAJOR REFACTOR - This file needs to be rewritten to remove all local_mode logic

---

#### 3. `src/tunacode/core/agents/agent_components/agent_config.py`

**Lines to modify:**
- **Line 25:** DELETE import of `is_local_mode` from `tunacode.utils.limits`
- **Lines 369-371:** DELETE local_mode check for minimal tool set
- **Lines 372-382:** DELETE minimal tool set definition
- **Lines 383-396:** KEEP - this becomes the standard tool set (remove the else wrapper)
- **Lines in `load_system_prompt()`:** DELETE local_mode template selection logic
- **Lines in agent creation:** DELETE local_mode `max_tokens` setting logic

**Action:** REMOVE conditional logic, keep full tool set as standard

---

#### 4. `src/tunacode/core/agents/resume/prune.py`

**Lines to modify:**
- **Line 11:** DELETE import of `is_local_mode`
- **Lines 17-18:** DELETE `LOCAL_PRUNE_PROTECT_TOKENS`, `LOCAL_PRUNE_MINIMUM_THRESHOLD`
- **Lines 30, 36-37:** DELETE `get_prune_thresholds()` function - replace with direct constants
- Update all callers to use standard thresholds directly

**Action:** DELETE function, use standard thresholds

---

#### 5. `src/tunacode/core/agents/resume/summary.py`

**Lines to modify:**
- **Line 28:** DELETE `LOCAL_SUMMARY_THRESHOLD` constant
- **Line 98:** DELETE `local_mode: bool = False` parameter
- **Line 104:** DELETE comment about local_mode parameter
- **Line 112:** DELETE conditional threshold selection

**Action:** Remove parameter, always use standard `SUMMARY_THRESHOLD`

---

#### 6. `src/tunacode/core/prompting/templates.py`

**Lines to modify:**
- **Lines 62-69:** DELETE `LOCAL_TEMPLATE` definition

**Action:** DELETE template definition

---

#### 7. `src/tunacode/core/prompting/__init__.py`

**Lines to modify:**
- **Line with `LOCAL_TEMPLATE` import:** DELETE from import

**Action:** Remove from `__all__`

---

#### 8. `src/tunacode/tools/bash.py`

**Lines to modify:**
- **Line 15:** DELETE import of `LOCAL_MAX_COMMAND_OUTPUT`
- **Lines 18-23:** DELETE `_get_max_output()` function
- **Line in truncation logic:** Use `MAX_COMMAND_OUTPUT` directly

**Action:** Remove local_mode conditional, use standard constant

---

#### 9. `src/tunacode/tools/read_file.py`

**Lines to modify:**
- **Lines 14-15:** DELETE imports of `LOCAL_DEFAULT_READ_LIMIT`, `LOCAL_MAX_LINE_LENGTH`
- **Lines 21-26:** DELETE `_get_limits()` function
- **Line in read logic:** Use `DEFAULT_READ_LIMIT` and `MAX_LINE_LENGTH` directly

**Action:** Remove local_mode conditional, use standard constants

---

### Test Files

#### `tests/unit/utils/test_limits.py`

**Status:** This entire test file is testing local_mode functionality.

**Lines affected:**
- Lines 3, 43, 49, 51, 55, 59, 73, 76, 89, 92, 105, 109, 130, 134, 155, 159, 181, 184, 189, 192, 197, 200, 223

**Action:** DELETE entire test file OR rewrite to test new simplified behavior

---

### Documentation Files

All documentation files reference `local_mode`. After removal, these should be updated:

1. **`docs/codebase-map/modules/INDEX.md`** - Update utils-limits description
2. **`docs/codebase-map/modules/core-agents.md`** - Remove local mode sections
3. **`docs/codebase-map/modules/core-compaction.md`** - Remove local mode sections
4. **`docs/codebase-map/modules/utils-limits.md`** - Complete rewrite or delete
5. **`docs/codebase-map/architecture/architecture.md`** - Remove local mode sections
6. **`docs/codebase-map/architecture/compaction-ontology.md`** - Remove local mode references
7. **`docs/configuration/README.md`** - Remove local mode documentation
8. **`docs/configuration/tunacode.json.example`** - Remove `local_mode` and `local_max_tokens` fields
9. **`docs/configuration/tunacode.local.json.example`** - DELETE this entire file

---

### Constants to Delete

```python
# From src/tunacode/constants.py
LOCAL_MAX_COMMAND_OUTPUT = 1500
LOCAL_DEFAULT_READ_LIMIT = 200
LOCAL_MAX_LINE_LENGTH = 500
LOCAL_MAX_FILES_IN_DIR = 20

# From src/tunacode/core/agents/resume/prune.py
LOCAL_PRUNE_PROTECT_TOKENS = 2000
LOCAL_PRUNE_MINIMUM_THRESHOLD = 500

# From src/tunacode/core/agents/resume/summary.py
LOCAL_SUMMARY_THRESHOLD = 6000
```

---

### Functions to Delete

```python
# From src/tunacode/utils/limits.py
def is_local_mode() -> bool:
    """Check if local_mode is enabled."""
    return _load_settings().get("local_mode", False)

# From src/tunacode/core/agents/resume/prune.py
def get_prune_thresholds() -> tuple[int, int]:
    """Get pruning thresholds based on local_mode setting."""
    # ... entire function
```

---

### Templates to Delete

```python
# From src/tunacode/core/prompting/templates.py
LOCAL_TEMPLATE = """{{AGENT_ROLE}}

====

{{TOOL_USE}}

====

{{USER_INSTRUCTIONS}}"""
```

---

## Key Patterns / Solutions Found

### Single Toggle, System-Wide Impact

The `local_mode` setting was a single boolean in config that cascaded through 6 optimization layers. This was implemented via:

1. **Central detection** (`utils/limits.py::is_local_mode()`)
2. **Import by consumers** (agent_config, prune, summary, bash, read_file)
3. **Conditional logic** in each module

### Replacement Path

Users can now achieve same goals by:
- Setting `OPENAI_BASE_URL` for local model endpoints
- Configuring individual limits explicitly in settings
- Using models with larger context windows

## Knowledge Gaps

1. **Current user base:** How many users rely on `local_mode`?
2. **Migration guide:** Should we provide a migration path?
3. **Default values:** After removal, what should the new defaults be?

## References

### GitHub Permalink to Original Commit
https://github.com/alchemiststudiosDOTai/tunacode/blob/cb5a4568a812696e7d41796736a498cf77b6ca96

### Production Code Files
- `src/tunacode/constants.py:34-38`
- `src/tunacode/utils/limits.py` (entire file)
- `src/tunacode/core/agents/agent_components/agent_config.py`
- `src/tunacode/core/agents/resume/prune.py`
- `src/tunacode/core/agents/resume/summary.py`
- `src/tunacode/core/prompting/templates.py:62-69`
- `src/tunacode/core/prompting/__init__.py`
- `src/tunacode/tools/bash.py`
- `src/tunacode/tools/read_file.py`

### Test Files
- `tests/unit/utils/test_limits.py` (entire file)
- `tests/test_compaction.py` (some tests)

### Documentation Files
- `docs/codebase-map/modules/INDEX.md`
- `docs/codebase-map/modules/core-agents.md`
- `docs/codebase-map/modules/core-compaction.md`
- `docs/codebase-map/modules/utils-limits.md`
- `docs/codebase-map/architecture/architecture.md`
- `docs/codebase-map/architecture/compaction-ontology.md`
- `docs/configuration/README.md`
- `docs/configuration/tunacode.json.example`
- `docs/configuration/tunacode.local.json.example` (DELETE)

### Related Cards
- None (this is pre-implementation research)
