---
name: 'Update 5 historical plan files still referencing source: .claude/BACKLOG.md'
description: '5 plan files in plan/ still have source: .claude/BACKLOG.md in their metadata. Files: tasks-2-validator-ux-coverage.md, architect-validate-orchestrator-discipline.md, feature-context-conventional-commits-changelog-refs.md, tasks-3-validator-ci-gate.md, tasks-4-validate-orchestrator-discipline.md. Update source references to point to the specific per-item file paths in .claude/backlog/ instead.'
metadata:
  topic: update-5-historical-plan-files-still-referencing-source-clau
  source: Workflow validation session 2026-02-27
  added: '2026-02-27'
  priority: P2
  type: Documentation
  status: in-progress
  issue: '#287'
  groomed: '2026-02-27'
  last_synced: '2026-02-27T12:04:55Z'
  plan: N/A
---

## Fact-Check

**Claims checked: 7 | VERIFIED: 3 | REFUTED: 4 | INCONCLUSIVE: 0**

| Claim | Verdict | Evidence |
|-------|---------|----------|
| 5 plan files have `source: .claude/BACKLOG.md` in metadata | REFUTED | grep found only 1 frontmatter `source:` + 1 markdown `**Source**:` referencing BACKLOG.md. Other 3 mention BACKLOG.md in content, not as source metadata |
| `tasks-2-validator-ux-coverage.md` has source reference | VERIFIED | Line 5: `source: .claude/BACKLOG.md (P1)` |
| `architect-validate-orchestrator-discipline.md` has source reference | REFUTED | No source field; BACKLOG.md appears only in test scenario table |
| `feature-context-conventional-commits-changelog-refs.md` has source reference | VERIFIED | Line 7: `**Source**: ...BACKLOG.md` |
| `tasks-3-validator-ci-gate.md` has source reference | REFUTED | No source metadata; BACKLOG.md in action items only |
| `tasks-4-validate-orchestrator-discipline.md` has source reference | REFUTED | No source metadata; BACKLOG.md in test scenarios only |
| `.claude/BACKLOG.md` is obsolete (replaced by per-item files) | VERIFIED | Backlog uses `.claude/backlog/` per-item files |

**Scope correction**: Actual work is 2 files (not 5). The other 3 files mention BACKLOG.md in body content discussing test scenarios — these are legitimate content references, not stale source pointers.

## RT-ICA

**Goal**: Remove stale references to the deleted .claude/BACKLOG.md in plan file metadata, replacing them with correct per-item backlog file paths.

| # | Condition | Status | Info |
|---|-----------|--------|------|
| 1 | List of affected files | AVAILABLE (scope corrected) | Only 2 of 5 claimed files have source metadata references: `tasks-2-validator-ux-coverage.md` (frontmatter), `feature-context-conventional-commits-changelog-refs.md` (markdown metadata) |
| 2 | Correct replacement paths | DERIVABLE | `tasks-2-validator-ux-coverage.md` → `.claude/backlog/p2-add-pr003pr004-test-coverage-to-plugin-registration-validato.md`; `feature-context-conventional-commits-changelog-refs.md` → `.claude/backlog/p2-conventional-commits-fix-changelog-references-to-nonexistent.md` |
| 3 | Whether .claude/BACKLOG.md still exists | DERIVABLE | Backlog system uses per-item files; verify with filesystem |
| 4 | Plan files not actively in use | DERIVABLE | Both show status complete/DISCOVERY_COMPLETE |

**Decision**: APPROVED
**Missing**: None — scope is 2 files (not 5) but all conditions met.

## Groomed (2026-02-27)

### Scope

**Corrected scope**: 2 files need source metadata updates (not 5 as originally claimed). The other 3 files mention BACKLOG.md in content body discussing test scenarios — these are legitimate content references, not stale source pointers.

### Priority

P2 — Technical debt cleanup post-backlog system migration. Minimal effort (2 lines across 2 files). Purely cosmetic/archival.

### Impact

- Corrects metadata sourcing to match current backlog architecture (per-item files in `.claude/backlog/`)
- Enables proper artifact tracing from plan files back to original backlog items
- Reduces confusion in future historical investigations

### Output / Evidence

1. `plan/tasks-2-validator-ux-coverage.md` line 5: update `source: .claude/BACKLOG.md (P1)` → `source: .claude/backlog/p2-add-pr003pr004-test-coverage-to-plugin-registration-validato.md`
2. `plan/feature-context-conventional-commits-changelog-refs.md` line 7: update `**Source**: Plugin code review session 2026-02-21, BACKLOG.md` → `**Source**: Plugin code review session 2026-02-21, .claude/backlog/p2-conventional-commits-fix-changelog-references-to-nonexistent.md`
3. Grep confirms no remaining `source:` or `**Source**:` fields referencing `.claude/BACKLOG.md` in plan/ metadata

### Dependencies

- Depends on: None — standalone cleanup
- Blocks: None — archival only

### Research

None required — all data verified via filesystem inspection.

### Skills

None required.

### Agents

None required — trivial 2-line edit.

### Prior Work

- Backlog system redesign (#282) migrated from monolithic BACKLOG.md to per-item files
- Both plan files are for completed work (status: complete / DISCOVERY_COMPLETE)

### Files

| File | Action |
|------|--------|
| `plan/tasks-2-validator-ux-coverage.md` | Edit line 5 frontmatter `source:` |
| `plan/feature-context-conventional-commits-changelog-refs.md` | Edit line 7 `**Source**:` |
| `.claude/backlog/p2-add-pr003pr004-test-coverage-to-plugin-registration-validato.md` | Target reference for file 1 |
| `.claude/backlog/p2-conventional-commits-fix-changelog-references-to-nonexistent.md` | Target reference for file 2 |

### Decision

APPROVED — Scope reduced from 5 files to 2. Effort: small (two text replacements). No architectural decisions needed.