# Technical Design: fix-task-list-id

## Metadata
- **Feature**: fix-task-list-id
- **Status**: DRAFT
- **Created**: 2026-01-31

---

## 1. Overview

### 1.1 Summary
Fix CLAUDE_CODE_TASK_LIST_ID propagation so solo-dev orchestrators and workers share the same default task list. Currently, workers get `feature` as fallback while orchestrators use the default list, causing task coordination to silently fail.

### 1.2 Goals
- Workers inherit orchestrator's CLAUDE_CODE_TASK_LIST_ID (or absence thereof)
- Teams can still opt-in by exporting the env var before starting Claude Code
- No behavioral change for team workflows that already set the env var

### 1.3 Non-Goals
- Changing the task list ID format or naming convention
- Adding new env vars or config options

---

## 2. Architecture

### 2.1 Root Cause

```
Orchestrator (no env var set)  →  uses DEFAULT task list
    │
    └── spawns worker with CLAUDE_CODE_TASK_LIST_ID=feature
            → uses FEATURE task list  ← MISMATCH
```

### 2.2 Fix Strategy

Remove `feature` fallback from all 4 spawn methods in launcher.py. Only propagate the env var if the orchestrator actually has it set. Remove the safety net in worker_main.py that forces it. Update rush.core.md container template.

### 2.3 After Fix

```
Solo:  orchestrator(default) → worker(default)     ✓
Team:  orchestrator(feat)    → worker(feat)         ✓
```

---

## 3. Detailed Changes

### 3.1 launcher.py — 4 spawn sites

Lines ~462-466, ~734-737 (subprocess):
```python
# Before:
worker_env.setdefault(
    "CLAUDE_CODE_TASK_LIST_ID",
    os.environ.get("CLAUDE_CODE_TASK_LIST_ID", feature),
)
# After:
task_list_id = os.environ.get("CLAUDE_CODE_TASK_LIST_ID")
if task_list_id:
    worker_env.setdefault("CLAUDE_CODE_TASK_LIST_ID", task_list_id)
```

Lines ~974-976, ~1543-1545 (container):
```python
# Before:
"CLAUDE_CODE_TASK_LIST_ID": os.environ.get("CLAUDE_CODE_TASK_LIST_ID", feature),
# After:
**({"CLAUDE_CODE_TASK_LIST_ID": os.environ["CLAUDE_CODE_TASK_LIST_ID"]}
   if "CLAUDE_CODE_TASK_LIST_ID" in os.environ else {}),
```

### 3.2 worker_main.py — line ~116-117
Delete:
```python
if not env.get("CLAUDE_CODE_TASK_LIST_ID"):
    env["CLAUDE_CODE_TASK_LIST_ID"] = args.feature
```

### 3.3 zerg:rush.core.md — line ~178
```bash
# Before:
-e CLAUDE_CODE_TASK_LIST_ID=$FEATURE \
# After:
-e CLAUDE_CODE_TASK_LIST_ID=${CLAUDE_CODE_TASK_LIST_ID:-} \
```

### 3.4 CLAUDE.md — doc updates
- Rule 5: note inheritance behavior
- Drift check #5: lower threshold from ≥5 to ≥4

---

## 4. Key Decisions

### Decision: Omit env var when unset (vs empty string)

**Context**: When orchestrator has no CLAUDE_CODE_TASK_LIST_ID, should workers get an empty string or no env var at all?

**Decision**: No env var (omit entirely)

**Rationale**: Empty string may be interpreted differently than absence by Claude Code internals. Omitting preserves the default behavior exactly.

---

## 5. Implementation Plan

| Phase | Tasks | Parallel |
|-------|-------|----------|
| Foundation (L1) | 2 tasks: launcher.py + worker_main.py | Yes |
| Integration (L2) | 1 task: rush.core.md + CLAUDE.md docs | No (depends on L1) |
| Verification (L3) | 1 task: run tests + grep validation | No (depends on L2) |

### File Ownership

| File | Task ID | Operation |
|------|---------|-----------|
| zerg/launcher.py | TASK-L1-001 | modify |
| zerg/worker_main.py | TASK-L1-002 | modify |
| zerg/data/commands/zerg:rush.core.md | TASK-L2-001 | modify |
| CLAUDE.md | TASK-L2-001 | modify |

---

## 6. Risk Assessment

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| Teams relying on implicit feature fallback | Low | Med | Doc update explains new behavior |
| Container env dict syntax change | Low | Low | Test with existing container tests |

---

## 7. Verification

```bash
# 1. All existing tests pass
python -m pytest tests/ -x -q

# 2. launcher.py no longer has feature fallback
grep -n "CLAUDE_CODE_TASK_LIST_ID" zerg/launcher.py
# Expect: ALLOWED_ENV_VARS entry + conditional sets (no `feature` fallback)

# 3. worker_main.py no longer forces the env var
grep -c "CLAUDE_CODE_TASK_LIST_ID" zerg/worker_main.py
# Expect: 0

# 4. Drift detection still passes (with updated threshold)
count=$(grep -c "CLAUDE_CODE_TASK_LIST_ID" zerg/launcher.py)
echo "launcher.py — $count CLAUDE_CODE_TASK_LIST_ID refs (expect ≥4)"
```
