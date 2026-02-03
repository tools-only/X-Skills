# fn-22-53k Epic-completion review gate

## Overview

Add review gate that runs when all epic tasks complete, BEFORE epic closes. Compares actual implementation against spec to find requirement gaps.

## Goal: Spec Compliance Only

**This review checks: "Does the implementation deliver everything the spec requires?"**

Verifies ALL deliverables mentioned in spec/acceptance criteria:
- Code changes implement required functionality
- Docs updated as specified (README, API docs, etc.)
- Tests added per acceptance criteria
- Config/schema changes if mentioned

Does NOT re-review code quality — that's handled by per-task impl-review.

**Catches:**
- Requirements that never became tasks (decomposition gaps)
- Requirements partially implemented across tasks (cross-task gaps)
- Scope drift (task marked done without fully addressing spec intent)
- Missing doc updates

**Verdicts:** `SHIP` / `NEEDS_WORK` (same as impl-review for consistency)

**Fix loop:** Same pattern as impl-review (works in both interactive and Ralph modes):
1. Review finds gaps → NEEDS_WORK
2. Agent implements missing requirements
3. Re-review (same chat for RP, same session for codex)
4. Repeat until SHIP
5. Receipt written, epic can close

In Ralph: fix loop happens inside skill invocation, not across iterations. Skill returns after SHIP.

**Gap task creation:** Only for complex gaps that truly need separate task scope. Primary behavior is fix-inline.

## Problem

Per-task impl-review validates task acceptance criteria, but:
- Task decomposition can miss requirements (spec items never become tasks)
- Task scope drifts (marked done without fully addressing intent)
- Cross-task requirements fall through cracks
- Manual review consistently finds 5-10 gaps per epic

## Solution

```
All tasks "done" → epic-completion-review → SHIP? close epic
                                         → NEEDS_WORK? fix gaps → re-review (loop until SHIP)
```

Two-phase approach (prevents over-correction bias per ASE'25 research):
1. Extract requirements from spec as explicit bullets (no code analysis yet)
2. Verify each requirement against actual code changes

## Naming Convention

**Consistent naming (mirrors plan_review pattern):**
- Epic field: `completion_review_status` (unknown | ship | needs_work)
- Epic field: `completion_reviewed_at` (timestamp)
- Selector status: `completion_review`
- Selector flag: `--require-completion-review`
- Ralph env: `COMPLETION_REVIEW=rp|codex|none` (controls whether review runs)
- Receipt type: `completion_review`
- Receipt filename: `completion-fn-N.json`

**Note:** No separate `REQUIRE_COMPLETION_REVIEW` env var — `COMPLETION_REVIEW != none` is the gate. Ralph passes `--require-completion-review` to selector when `COMPLETION_REVIEW != none`.

## Scope

**In scope:**
- `flowctl codex completion-review fn-N` command (LLM-driven review)
- `flowctl next --require-completion-review` flag
- `/flow-next:epic-review` skill (rp/codex backends)
- Command shim: `commands/flow-next/epic-review.md`
- Ralph gate in `maybe_close_epics()` (receipt + status check)
- Interactive `/flow-next:work` handling of `completion_review` status

**Out of scope:**
- Changing impl-review/plan-review behavior
- Standalone `flowctl epic-review` (verdict requires LLM; use skill or codex command)

## Quick commands

```bash
./plugins/flow-next/scripts/smoke_test.sh
./plugins/flow-next/scripts/ralph_smoke_test.sh
```

## Acceptance

- [ ] `flowctl next --require-completion-review` returns `status=completion_review` when epic needs review
- [ ] `flowctl codex completion-review fn-N` runs LLM review and outputs verdict
- [ ] Ralph gates closure on BOTH `completion_review_status=ship` AND valid receipt
- [ ] `/flow-next:epic-review` skill works with rp/codex backends
- [ ] Command shim `epic-review.md` invokes skill
- [ ] Interactive `/flow-next:work` handles `completion_review` status
- [ ] Works in both Ralph and interactive modes

## Integration Points

| Component | File | Change |
|-----------|------|--------|
| Epic JSON | `flowctl.py:598` | Add `completion_review_status`, `completion_reviewed_at` in `normalize_epic()` |
| Selector | `flowctl.py:4072` | Add `--require-completion-review` flag, return `completion_review` status |
| CLI | `flowctl.py` | New `epic set-completion-review-status` command |
| CLI | `flowctl.py` | New `codex completion-review` command |
| Ralph | `ralph.sh:751` | Gate in `maybe_close_epics()`: check status AND verify receipt |
| Ralph | `ralph.sh:876` | Add `elif [[ "$status" == "completion_review" ]]` dispatch |
| Ralph | `ralph.sh` | New `prompt_completion.md` template |
| Hook | `ralph-guard.py` | Handle `completion-fn-N.json` receipt pattern |
| Skill | new | `flow-next-epic-review/SKILL.md` + `workflow.md` |
| Command | new | `commands/flow-next/epic-review.md` (invokes skill) |
| Work skill | existing | Update `flow-next-work/phases.md` for `completion_review` status |

## References

- Issue: https://github.com/gmickel/gmickel-claude-marketplace/issues/83
- impl-review pattern: `plugins/flow-next/skills/flow-next-impl-review/`
- plan-review pattern: `plugins/flow-next/skills/flow-next-plan-review/`
- maybe_close_epics: `ralph.sh:751-771`
- selector: `flowctl.py:4000-4171`
