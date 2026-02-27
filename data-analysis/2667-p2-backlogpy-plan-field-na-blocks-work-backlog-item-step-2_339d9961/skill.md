---
name: 'backlog.py: plan field N/A blocks work-backlog-item Step 2'
description: "When per-item files have plan: N/A (set on completed items), work-backlog-item Step 2 treats it as a valid plan path and stops with 'This item already has a plan at N/A'. Currently all 22 items with plan: N/A are status: done so they're filtered out, but this is a latent bug. Fix: either strip N/A values during parsing in backlog.py, or add a guard in Step 2 to ignore non-path plan values."
metadata:
  topic: backlogpy-plan-field-na-blocks-work-backlog-item-step-2
  source: Workflow validation session 2026-02-27
  added: '2026-02-27'
  priority: P2
  type: Bug
  status: open
  issue: '#285'
  groomed: '2026-02-27'
  last_synced: '2026-02-27T12:06:05Z'
---

## Fact-Check

**Date**: 2026-02-27
**Claims checked**: 4

| # | Claim | Verdict | Evidence |
|---|-------|---------|----------|
| 1 | Per-item files have plan: N/A set on completed items | VERIFIED | 22 backlog files have plan: N/A in frontmatter, all with status: done |
| 2 | Step 2 treats N/A as a valid plan path and stops | VERIFIED | work-backlog-item/SKILL.md:199-205 checks plan field truthiness; "N/A" is truthy |
| 3 | All items with plan: N/A are status: done | VERIFIED | All 22 files confirmed status: done, priority: completed |
| 4 | Description says 22 items | VERIFIED | Count is 22 (23rd grep match is bug report description text, not frontmatter) |

**Summary**: VERIFIED: 4 | REFUTED: 0 | INCONCLUSIVE: 0

## RT-ICA

**Goal**: Prevent work-backlog-item Step 2 from treating "N/A" plan values as valid plan paths.

| # | Condition | Status | Info |
|---|-----------|--------|------|
| 1 | Bug location identified | AVAILABLE | backlog.py:260 builds plan field; work-backlog-item/SKILL.md:199-205 checks it |
| 2 | Root cause understood | AVAILABLE | str("N/A") is truthy, treated as valid path |
| 3 | Fix strategy clear | AVAILABLE | Strip N/A during parsing OR add guard in Step 2; both viable |
| 4 | Files affected by fix identified | DERIVABLE | backlog.py (parsing), work-backlog-item/SKILL.md (Step 2 logic), or both |
| 5 | Regression scope known | AVAILABLE | 22 items with plan: N/A, all status:done |
| 6 | Test approach defined | DERIVABLE | Verify backlog list JSON output for plan field handling |
| 7 | Interaction with close/update/sync | DERIVABLE | Check if close also sets plan: N/A |

**Decision**: APPROVED
**Missing**: None — all conditions AVAILABLE or DERIVABLE

## Groomed (2026-02-27)

### Reproducibility

1. Run `uv run .claude/skills/backlog/scripts/backlog.py list --format json`
2. Filter results for items with `status: done` or `status: resolved`
3. Note items with `plan: "N/A"` in JSON output
4. Invoke `/work-backlog-item` with title of one of these items
5. Observe Step 2 checks `if plan field is truthy` (SKILL.md line 199-205)
6. Confirm false-positive: "This item already has a plan at N/A"

### Output / Evidence

- `backlog.py` line 260: `item["**Plan**"] = str(meta.get("plan") or fm.get("plan") or "")`
  - Converts `plan: N/A` to string `"N/A"`, which is truthy in Python
- `backlog.py` line 1003: `_update_item_metadata(..., {"metadata": {"status": "done", ... "plan": plan}})`
  - Explicitly sets `plan: N/A` when closing items
- `work-backlog-item/SKILL.md` lines 199-205: checks `if plan field is truthy` without filtering out `"N/A"`
- Current barrier: 22 items with `plan: N/A` all have `status: done`, so they're filtered at JSON stage; latent bug surfaces if a completed item gets reopened or a non-completed item receives `plan: N/A`

### Priority

7/10 — Prevents `work-backlog-item` from processing items marked done (current safety by accident), but exposes code smell: N/A should normalize to empty, not stringify to truthy value. Low immediate risk (all affected items filtered), moderate debt (confuses workflow state).

### Impact

- Blocks: Items with `plan: N/A` cannot be reopened and re-planned without manual frontmatter edit
- Bottleneck: Design confusion between "completed" (status: done) and "no plan needed" (plan: N/A); both conflate in Step 2 logic

### Scope

Two viable fixes (both produce same outcome):

**Option A: Normalize in backlog.py parsing (recommended)**
- Line 260: Change `str(meta.get("plan") or fm.get("plan") or "")` to filter `N/A` to empty
- Applies globally; all downstream code sees clean data
- Single source of truth

**Option B: Guard in work-backlog-item Step 2**
- Lines 199-205: Add condition `if plan and plan != "N/A"` before stopping
- Localized; other skills/agents using backlog.py still receive N/A strings

### Output / Evidence

**Expected Behavior**: When `work-backlog-item` Step 2 extracts plan field, it should treat `"N/A"` as "no plan available" (equivalent to empty string), not as a valid plan path.

**Acceptance Criteria**:
1. `backlog.py list --format json` on items with `plan: N/A` returns `plan: ""` (empty string)
2. `/work-backlog-item` with a reopened item does NOT report "plan at N/A"
3. Existing 22 files with `plan: N/A` are unaffected (status: done)
4. New items closing without a plan serialize to `plan: ""` in JSON

### Dependencies

- Depends on: None
- Blocks: None (safe to do independently)

### Research

No external research needed. Bug is fully characterized from codebase analysis.

### Skills

- `/work-backlog-item` — consumer of the plan field
- `/backlog` — manages per-item files and JSON serialization

### Agents

None required for fix implementation.

### Prior Work

- `backlog.py` lines 260, 1003 — plan field serialization and close logic
- `work-backlog-item/SKILL.md` lines 199-205 — Step 2 plan guard

### Files

- `.claude/skills/backlog/scripts/backlog.py` (lines 260, 1003)
- `.claude/skills/work-backlog-item/SKILL.md` (lines 199-205)

### Decision

APPROVED — Small fix, no blockers. Recommend Option A (normalize at parse time).

### Effort

Small — 1-line string filter or 1-line guard condition. ~15 minutes implementation + validation.