---
feature: validator-ci-gate
slug: validator-ci-gate
created: 2026-02-18
status: complete
backlog-item: "Resolve plugin-validator pre-existing errors to make CI gate blocking"
priority: P1
---

# Task Plan: Make validate-plugins CI gate blocking

## Background

The backlog item (added 2026-02-14) described 19 errors and 63 warnings preventing
the validate-plugins CI step from being blocking. As of 2026-02-18, the situation
has changed dramatically:

**Current validator state (verified 2026-02-18)**:
- Exit code: 0 (no errors)
- 0 LK001 errors (broken links — all fixed in prior work)
- 0 SK007 errors (orchestrating-swarms already split into 4 sub-skills)
- 0 SK005 warnings (trigger phrases — resolved)
- 13 SK006 warnings (complexity — these are warnings, not errors)
- 4 LK002 warnings (missing `./` prefix in agentskills/SKILL.md)

The validator already passes. The one-line CI change to make it blocking is safe now.

## Task 1: Fix LK002 warnings in agentskills/SKILL.md

- **Status**: COMPLETE
- **Dependencies**: none
- **Priority**: low
- **Complexity**: trivial
- **Agent**: general-purpose
- **Completed**: 2026-02-18

Add `./` prefix to 4 bare relative links in
`plugins/plugin-creator/skills/agentskills/SKILL.md`:
- `references/best-practices.md` → `./references/best-practices.md`
- `references/specification.md` → `./references/specification.md`
- `references/integration.md` → `./references/integration.md`

(Note: `best-practices.md` appears twice in the file.)

### Acceptance Criteria

1. All 4 LK002 warnings eliminated for agentskills/SKILL.md
2. Links still resolve to correct files
3. No new validator findings introduced

### Verification Steps

1. Run validator on the single file: `uv run plugins/plugin-creator/scripts/plugin_validator.py plugins/plugin-creator/skills/agentskills/SKILL.md`
2. Confirm 0 LK002 warnings
3. Verify referenced files exist at the `./references/` paths

## Task 2: Promote validate-plugins to blocking CI gate

- **Status**: COMPLETE
- **Dependencies**: none (validator already exits 0)
- **Priority**: high
- **Complexity**: trivial
- **Agent**: general-purpose
- **Completed**: 2026-02-18

Remove `validate-plugins` from the `allowed-failures` list at line 493 of
`.github/workflows/code-quality.yml`. Also update the comment block at lines
207-212 to reflect the new blocking status.

### Acceptance Criteria

1. `validate-plugins` is NOT in `allowed-failures`
2. Quality gate `needs` array still includes `validate-plugins`
3. Comment block accurately reflects blocking status
4. YAML syntax valid

### Verification Steps

1. Validate YAML: `python3 -c "import yaml; yaml.safe_load(open('.github/workflows/code-quality.yml'))"`
2. Run `uv run prek run --files .github/workflows/code-quality.yml`
3. After push: `gh run list -R Jamie-BitFlight/claude_skills --limit=3` to verify CI passes

## Task 3: Update backlog item as completed

- **Status**: COMPLETE
- **Dependencies**: Task 1, Task 2
- **Priority**: medium
- **Complexity**: trivial
- **Agent**: orchestrator (self)
- **Completed**: 2026-02-18

Move the backlog item from P1 to Completed section in `.claude/BACKLOG.md`.
Update `p1-count` in frontmatter. Add completion date and description of what
was done.

### Acceptance Criteria

1. Item moved to Completed section with date and summary
2. Frontmatter counts accurate
3. No formatting issues in BACKLOG.md

### Verification Steps

1. Grep for the item title in Completed section
2. Verify p1-count decremented
3. Validate markdown structure

## Out of Scope

- **SK006 warnings (13 skills)**: These are complexity warnings, not errors.
  They don't block the CI gate. Each is a separate refactoring effort tracked
  as future work. The backlog item for "plugin-validator pre-commit output is
  too noisy" (P1) covers the UX concern around warning volume.
- **SK005 trigger phrase improvements**: Already resolved. If re-introduced,
  tracked separately.
- **agent-browser missing references**: The skill has only SKILL.md with no
  references/ or templates/ dirs. The broken links were apparently cleaned up
  in a prior commit. If re-introduced, tracked separately.

## Notes

The investigation report at
`plugins/plugin-creator/.claude/reports/validator-findings-investigation.md`
was comprehensive but is now partially stale — many of its findings were
resolved between 2026-02-14 and 2026-02-18. The current validator run is the
authoritative source of truth.
