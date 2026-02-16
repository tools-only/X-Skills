# Technical Design: Context Engineering Structural Guardrails

## Metadata
- **Feature**: context-engineering-guardrails
- **Issue**: #39
- **Status**: DRAFT
- **Created**: 2026-02-01

---

## 1. Overview

### 1.1 Summary
Create automated enforcement of context engineering principles and Task ecosystem integrity rules. A single validation module codifies the 5 drift-detection checks from CLAUDE.md, adds split-pair consistency checks, and auto-splits oversized commands. Integrated into CI and pre-commit.

### 1.2 Goals
- Automate CLAUDE.md drift detection checklist as runnable validation
- Auto-split command files that exceed 300-line threshold via existing `CommandSplitter`
- Provide a command template so new commands inherit CE patterns by default
- Fail CI on violations; run in pre-commit for local enforcement

### 1.3 Non-Goals
- Token budget counting per command (fragile, changes with content)
- New plugin ABCs or validator framework classes
- Bracketed subject prefix format checking (convention, not structural)

---

## 2. Architecture

### 2.1 New Module: `zerg/validate_commands.py`

Single module, no new classes. Functions return `(bool, list[str])` matching `validation.py` pattern.

**5 validation checks:**

| Function | What it does | Source |
|----------|-------------|--------|
| `validate_task_references(commands_dir)` | Every base `.md` has TaskCreate/TaskUpdate/TaskList/TaskGet | CLAUDE.md drift #1 |
| `validate_backbone_depth(commands_dir)` | worker/status/merge/stop/retry have >=3 Task refs | CLAUDE.md drift #2 |
| `validate_split_pairs(commands_dir)` | .core.md <-> .details.md <-> parent .md all exist | CE integrity |
| `validate_split_threshold(commands_dir, auto_split)` | Files >=300 lines without .core.md get auto-split | CE enforcement |
| `validate_state_json_without_tasks(commands_dir)` | Files referencing .zerg/state must also ref TaskList/TaskGet | CLAUDE.md drift #3 |

Plus:
- `validate_all(commands_dir, auto_split)` — aggregator, runs all 5
- `_get_base_command_files(commands_dir)` — helper, returns .md files excluding .core.md/.details.md/_template.md
- `__main__` block for `python -m zerg.validate_commands`

**Constants:**
```python
BACKBONE_COMMANDS = {"worker", "status", "merge", "stop", "retry"}
BACKBONE_MIN_REFS = 3
TASK_MARKERS = {"TaskCreate", "TaskUpdate", "TaskList", "TaskGet"}
STATE_PATTERNS = re.compile(r"state.*json|STATE_FILE|\.zerg/state")
EXCLUDED_PREFIXES = ("_",)  # Skip _template.md
```

Imports `MIN_LINES_TO_SPLIT` from `command_splitter.py` (single source of truth).

### 2.2 Auto-Split Behavior

`validate_split_threshold` accepts `auto_split: bool = False`:
- When `False` (CI mode): returns error listing oversized files
- When `True` (pre-commit mode): calls `CommandSplitter.split_file()` on each oversized file, reports what was split, returns success

CI workflow detects uncommitted changes after validation and fails with message: "Command files were auto-split. Commit the generated files."

### 2.3 Command Template: `zerg/data/commands/_template.md`

~45 lines. Scaffold for new commands with:
- `# ZERG {CommandName}` header placeholder
- Arguments/flags section
- Pre-flight checks section
- Task Tracking boilerplate (pattern from `build.md:50-65`)
- Execution section placeholder
- Exit Codes section
- Help section

Underscore prefix auto-excludes from validator.

### 2.4 CI Integration: `.github/workflows/command-validation.yml`

New workflow (separate from changelog-check.yml):
- Triggers: push to main, PR events
- Steps: checkout, setup-python, pip install -e ., python -m zerg.validate_commands
- Fails build on any violation

### 2.5 Pre-commit Integration

Append local hook to `.pre-commit-config.yaml`:
```yaml
- repo: local
  hooks:
    - id: validate-commands
      name: Validate ZERG command files
      entry: python -m zerg.validate_commands --auto-split
      language: system
      files: 'zerg/data/commands/.*\.md$'
      pass_filenames: false
```

---

## 3. File Manifest

| File | Action | Est. Lines |
|------|--------|-----------|
| `zerg/validate_commands.py` | CREATE | ~180 |
| `tests/unit/test_validate_commands.py` | CREATE | ~220 |
| `zerg/data/commands/_template.md` | CREATE | ~45 |
| `.github/workflows/command-validation.yml` | CREATE | ~30 |
| `.pre-commit-config.yaml` | MODIFY | +8 |
| `docs/context-engineering.md` | MODIFY | +20 |
| `CLAUDE.md` | MODIFY | +5 |
| `CHANGELOG.md` | MODIFY | +2 |

---

## 4. Test Strategy

`tests/unit/test_validate_commands.py` — ~20 tests across 6 classes.

**Positive tests** (real `zerg/data/commands/` directory):
- All commands pass task reference check
- All backbone commands pass depth check
- All split pairs are consistent
- No state-json-without-tasks violations
- `validate_all()` passes clean on current codebase

**Negative tests** (tmp_path fixtures):
- Missing task reference flagged
- Shallow backbone file flagged
- Orphaned .core.md flagged
- Orphaned .details.md flagged
- Oversized unsplit file flagged
- State JSON without TaskList flagged
- _template.md excluded from checks
- .core.md/.details.md excluded from base file checks

**Auto-split tests:**
- Oversized file gets split when auto_split=True
- Split produces .core.md and .details.md
- Validator returns success after auto-split

---

## 5. Verification

```bash
# Validator passes on current codebase
python -m zerg.validate_commands

# New tests pass
pytest tests/unit/test_validate_commands.py -v

# Full suite regression
pytest tests/ --no-header -q

# Pre-commit hook runs
pre-commit run validate-commands --all-files
```
