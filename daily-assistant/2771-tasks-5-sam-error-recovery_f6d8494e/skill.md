---
description: "Research and design SAM error recovery / rollback procedures — cross-framework pattern synthesis, tier-based procedure design, and draft Appendix for stateless-software-engineering-framework.md"
version: "1.0"
tasks:
  - T1: Synthesize cross-framework error recovery patterns
  - T2: Design SAM error recovery procedure (3-tier model)
  - T3: Draft Appendix content for stateless-software-engineering-framework.md
  - T4: Integrate Appendix into external repo (depends on T1 + T2 + T3)
task_exports:
  enabled: false
  directory: "TASK"
---

# Task Plan: SAM Error Recovery / Rollback Procedures

## Context

**Backlog item**: SAM: Error Recovery / Rollback Procedures
**Priority**: P1
**Added**: 2026-02-01
**Source**: Gap analysis of SAM framework
**RT-ICA**: APPROVED (grooming-2026-02-21.md — 4 conditions, all DERIVABLE from research)
**Target**: New Appendix or Part 6 addition to [`stateless-software-engineering-framework.md`](https://github.com/bitflight-devops/stateless-agent-methodology/blob/main/stateless-software-engineering-framework.md)

**Problem statement**: SAM currently has no explicit procedure for when a task fails irrecoverably. There is no defined way to undo artifact changes or restore the artifact plane to a consistent state after a partial failure.

**Research questions** (from backlog):
1. How do GSD, BMAD-METHOD, AutoGPT, and traditional CI/CD handle rollback?
2. What patterns exist for transactional artifact updates?

---

## Context Manifest

### Research Findings (T1 — pre-synthesized from local research files)

The following findings are derived from research files available in this repository. They inform T2 (procedure design) and should be cited in T3 (Appendix draft).

#### GSD: Atomic Git Commits as Rollback Primitives

Source: `research/agent-frameworks/get-shit-done.md`

GSD implements one atomic git commit per completed task immediately upon task completion. This provides:
- **Granular rollback points**: `git bisect` locates failures at task granularity
- **Independent revertability**: each task commit can be reverted without affecting preceding or succeeding tasks
- **Session continuity**: future Claude sessions read clear, labeled history

The `gsd-debugger` agent is specialized for systematic debugging with state persistence. It maintains `STATE.md` across sessions so debug context survives context compaction.

**Import for SAM**: Adopt the task-level atomic commit discipline as the prerequisite for rollback. SAM agents should commit after each task completion — not after phase completion. This converts irrecoverable failures into recoverable ones: if Task N fails, revert the Task N commit.

**Gap in GSD for SAM**: GSD's rollback is always by `git revert` — it does not address in-memory artifact plane restoration (e.g., when an agent has produced malformed content in a SKILL.md before any commit). SAM needs a pre-commit artifact integrity check.

#### ARL/Ralph: Resource Bounds and Escalation Thresholds

Source: `plugins/plugin-creator/skills/arl/references/synthesis-arl-applicable.md` (R2 — Loop Detection)

Ralph implements two specific patterns relevant to SAM error recovery:
- **Count-based failure escalation** (`mod.rs:447, :1917`): After N task-level failures, escalate rather than continue. Ralph's default is abandon after 3 failures.
- **MaxRuntime cutoff** (`mod.rs:429-439`): Time-based termination as a safety net independent of failure count.

The ARL synthesis notes: "Task-based thrashing detection — Ralph counts task-level failures; ARL needs finding-level pattern detection." For SAM, the equivalent is phase-level failure detection.

**Import for SAM**: Define explicit escalation thresholds: N consecutive task failures → escalate to human. M total failures in a phase → abort phase and rollback all phase artifacts.

#### ARL General Theory: Structure Over Instruction

Source: `plugins/plugin-creator/skills/arl/references/synthesis-general-theory.md` (Principle 1)

"Telling an AI agent 'please do X' is unreliable. Structuring the pipeline so X is the only possible path is reliable." (SAM: ssf:108 — "Behavioral instructions cannot override architectural limitations.")

**Import for SAM**: Error recovery procedures must be structural (built into the pipeline), not instructional (agent told to recover). A RECOVERY gate in the SAM pipeline is reliable; an instruction "if something fails, clean up" is not.

#### CI/CD Patterns: Transactional Artifact Updates

Traditional CI/CD handles rollback through three complementary mechanisms:
1. **Idempotent operations**: Each task can be re-run without side effects (CREATE-OR-REPLACE, not APPEND-THEN-MODIFY).
2. **Staged commits**: Changes are staged in a branch/workspace, then applied atomically only on success. Partial completion leaves no artifacts.
3. **Health checks before promotion**: Before merging to main, verify artifact integrity. If health check fails, discard the branch workspace.

**Import for SAM**: The artifact plane (task files, SKILL.md files, plugin.json) should be treated as a staged workspace. Phase execution writes to a staging area. Only on phase SUCCESS are changes committed. On phase FAILURE, the staging area is discarded.

---

### Key Decisions (Do Not Re-Investigate)

**Decision 1 — Three-tier failure model**: SAM errors fall into three tiers based on recoverability:

| Tier | Name | Description | Recovery action |
|------|------|-------------|-----------------|
| T-1 | Recoverable | Task output is wrong but artifacts are intact; no commit yet | Retry the task |
| T-2 | Partial | Some artifacts modified, phase incomplete | Revert modified artifacts to pre-phase state, escalate |
| T-3 | Irrecoverable | Cannot determine artifact state, or commit already pushed | Manual review required; document in STATE.md |

**Decision 2 — Pre-commit integrity gate**: Before each task commit, run a lightweight structural check on modified artifacts (YAML frontmatter valid, required fields present, file not empty). If the gate fails, do not commit — apply T-1 recovery (retry).

**Decision 3 — Phase-level rollback boundary**: The rollback boundary is the phase, not the task. If a phase fails, revert ALL artifact changes from that phase using `git stash` or `git reset --soft HEAD~N` where N is the number of commits made in the current phase. The pre-phase state becomes the recovery baseline.

**Decision 4 — STATE.md failure record**: Every T-2 or T-3 failure is recorded in STATE.md with:
- Timestamp
- Phase and task where failure occurred
- Artifacts modified (list of file paths)
- Failure reason
- Recovery action taken

This provides the audit trail needed for human review and for future sessions to avoid repeating the same failure.

**Decision 5 — Appendix placement**: The procedure goes in a new `## Appendix B: Error Recovery and Rollback Procedures` section in `stateless-software-engineering-framework.md`, after the existing Appendix A (if any) or after the final Phase section. It is self-contained and cross-referenced from the Phase execution sections (3.x).

---

## T1: Synthesize Cross-Framework Error Recovery Patterns

**Status**: ✅ PRE-DONE (research synthesized above in Context Manifest)
**Dependencies**: None
**Priority**: 1
**Complexity**: Medium
**Agent**: Research synthesis (general-purpose)

**Target**: Research synthesis documented in Context Manifest above + local reference file
**Issue Type**: RESEARCH

**Description**:

Synthesize error recovery and rollback patterns from GSD, BMAD, Ralph/ARL, and CI/CD frameworks. Key patterns already synthesized above. This task produces a reference file that T3 cites in the Appendix draft.

**Acceptance Criteria**:
1. Local reference file `plan/sam-error-recovery-research.md` created with cross-framework comparison table
2. At least 4 frameworks covered: GSD, BMAD, Ralph/ARL, CI/CD
3. Each pattern includes: source, mechanism, applicability to SAM
4. Gaps (patterns missing from SAM) explicitly identified

**Required Inputs**:
- `research/agent-frameworks/get-shit-done.md`
- `research/agent-frameworks/bmad-method.md`
- `plugins/plugin-creator/skills/arl/references/synthesis-arl-applicable.md` (R2, R7)
- `plugins/plugin-creator/skills/arl/references/synthesis-general-theory.md` (Principles 1, 5)

**Expected Outputs**:
- `plan/sam-error-recovery-research.md` — cross-framework synthesis reference

**Can Parallelize With**: None (T2 depends on T1)

**Verification Steps**:
1. `wc -l plan/sam-error-recovery-research.md` — should be >50 lines
2. Verify 4+ frameworks mentioned: `grep -c "GSD\|BMAD\|Ralph\|CI/CD\|Temporal\|Airflow" plan/sam-error-recovery-research.md`
3. Verify gaps section exists: `grep "Gap\|Missing\|SAM lacks" plan/sam-error-recovery-research.md`

---

## T2: Design SAM Error Recovery Procedure

**Status**: ✅ PRE-DONE (decisions recorded above in Context Manifest)
**Dependencies**: T1
**Priority**: 1
**Complexity**: High
**Agent**: Architecture design (general-purpose)

**Target**: `plan/sam-error-recovery-design.md`
**Issue Type**: DESIGN

**Description**:

Design the error recovery procedure for SAM using the 3-tier failure model (T-1 Recoverable, T-2 Partial, T-3 Irrecoverable). The design must be:
- Structural (not instructional) per General Theory Principle 1
- Phase-scoped (rollback boundary = phase, not task)
- State-preserving (every failure recorded in STATE.md)
- Pre-commit gated (integrity check before each task commit)

Key design decisions already made (see Context Manifest). Do not reopen them.

**Acceptance Criteria**:
1. `plan/sam-error-recovery-design.md` documents the 3-tier failure model with decision criteria for each tier
2. Pre-commit integrity gate specification: what fields to check, which file types, exit conditions
3. Phase rollback procedure: exact git commands, scope definition, artifact list construction
4. STATE.md failure record format: all required fields, example entry
5. Cross-references to specific SAM SSF sections where each procedure integrates (e.g., Phase 3 execution, Phase 7 verification)

**Required Inputs**:
- `plan/sam-error-recovery-research.md` (T1 output)
- Cross-framework findings in Context Manifest
- 3-tier failure model (Decision 1 above)
- Phase-level rollback boundary (Decision 3 above)

**Expected Outputs**:
- `plan/sam-error-recovery-design.md` — detailed procedure specification

**Can Parallelize With**: None (T3 depends on T2)

**Verification Steps**:
1. Verify 3 tiers documented: `grep -c "T-1\|T-2\|T-3\|Recoverable\|Partial\|Irrecoverable" plan/sam-error-recovery-design.md`
2. Verify git commands present: `grep "git stash\|git reset\|git revert" plan/sam-error-recovery-design.md`
3. Verify STATE.md format defined: `grep "STATE.md" plan/sam-error-recovery-design.md`

---

## T3: Draft Appendix for stateless-software-engineering-framework.md

**Status**: ❌ NOT STARTED
**Dependencies**: T1 + T2
**Priority**: 1
**Complexity**: High
**Agent**: Documentation writer (general-purpose)

**Target**: `plan/sam-error-recovery-appendix-draft.md` (local staging file before cross-repo integration)
**Issue Type**: DOCUMENTATION

**Description**:

Write the full content for `## Appendix B: Error Recovery and Rollback Procedures` to be added to `stateless-software-engineering-framework.md` in `bitflight-devops/stateless-agent-methodology`.

The Appendix must:
1. Define the 3-tier failure classification model with decision tree
2. Specify the pre-commit integrity gate procedure
3. Specify the phase-level rollback procedure (git commands, artifact scope)
4. Define the STATE.md failure record format
5. Provide worked examples for each failure tier
6. Include cross-references to Phase sections in the SSF that invoke these procedures

Write in the existing SSF document style: structured sections, code blocks for git commands and file formats, tables for decision criteria.

**Acceptance Criteria**:
1. `plan/sam-error-recovery-appendix-draft.md` contains complete, standalone Appendix content
2. 3-tier model documented with at least one worked example per tier
3. Pre-commit integrity check procedure includes: which fields to validate, sample check script/pseudocode, pass/fail criteria
4. Phase rollback procedure includes: pre-conditions, exact git commands with flags, post-conditions (what state the artifact plane should be in), verification step
5. STATE.md failure record schema defined with all required fields and an example entry
6. Cross-reference list identifies ≥3 specific SSF sections that need to add `→ See Appendix B` links
7. No fabricated references to SSF section numbers — use `[TODO: confirm section N.N]` placeholders where the exact section is unknown
8. Citations for all cross-framework pattern claims (GSD research file, ARL synthesis files, CI/CD references)

**Required Inputs**:
- `plan/sam-error-recovery-research.md` (T1 output)
- `plan/sam-error-recovery-design.md` (T2 output)
- Existing SSF document structure: clone `bitflight-devops/stateless-agent-methodology` and read `stateless-software-engineering-framework.md` to match style

**Expected Outputs**:
- `plan/sam-error-recovery-appendix-draft.md` — complete Appendix content ready for cross-repo integration

**Can Parallelize With**: None (T4 depends on T3)

**Verification Steps**:
1. Verify draft exists and has content: `wc -l plan/sam-error-recovery-appendix-draft.md` — should be >100 lines
2. Verify 3 tiers documented: `grep -c "T-1\|T-2\|T-3\|Tier 1\|Tier 2\|Tier 3" plan/sam-error-recovery-appendix-draft.md`
3. Verify no fabricated SSF section numbers: `grep -v "TODO.*confirm" plan/sam-error-recovery-appendix-draft.md | grep "section [0-9]\." | head -5` — any results here need verification against actual SSF
4. Verify git commands are present: `grep "git\b" plan/sam-error-recovery-appendix-draft.md`
5. Verify citations block present: `grep "Source\|Citation\|research/" plan/sam-error-recovery-appendix-draft.md`

---

## T4: Integrate Appendix into External Repo

**Status**: ❌ NOT STARTED
**Dependencies**: T1 + T2 + T3
**Priority**: 2
**Complexity**: Medium
**Agent**: Cross-repo contributor (general-purpose)

**Target**: `stateless-software-engineering-framework.md` in `bitflight-devops/stateless-agent-methodology`
**Issue Type**: CROSS_REPO_INTEGRATION

**Description**:

1. Clone `https://github.com/bitflight-devops/stateless-agent-methodology`
2. Read existing `stateless-software-engineering-framework.md` to confirm section structure and find the correct insertion point (after last Part/Appendix section)
3. Replace all `[TODO: confirm section N.N]` placeholders in `plan/sam-error-recovery-appendix-draft.md` with actual section references from the SSF
4. Add cross-reference lines `→ See Appendix B: Error Recovery and Rollback Procedures` to the Phase execution sections identified in T3
5. Insert the complete Appendix content at the correct position
6. Create a PR against the external repo

**Acceptance Criteria**:
1. `stateless-software-engineering-framework.md` in external repo has the new Appendix B section
2. All `[TODO]` placeholders replaced with actual section references
3. At least 3 inbound cross-references from Phase execution sections to Appendix B
4. PR created with conventional-commits subject: `docs(ssf): add Appendix B — Error Recovery and Rollback Procedures`
5. PR body links back to this plan file and the originating backlog item

**Required Inputs**:
- `plan/sam-error-recovery-appendix-draft.md` (T3 output)
- Access to `bitflight-devops/stateless-agent-methodology` (requires GitHub write permissions)

**Expected Outputs**:
- PR in `bitflight-devops/stateless-agent-methodology` with Appendix B integration

**Pre-condition check before starting**:
```bash
gh repo clone bitflight-devops/stateless-agent-methodology /tmp/sam-repo
ls /tmp/sam-repo/stateless-software-engineering-framework.md
wc -l /tmp/sam-repo/stateless-software-engineering-framework.md
```

**Verification Steps**:
1. Appendix B heading present: `grep "Appendix B" /tmp/sam-repo/stateless-software-engineering-framework.md`
2. No TODO placeholders remain: `grep "TODO.*confirm" /tmp/sam-repo/stateless-software-engineering-framework.md` — should return nothing
3. Cross-references present: `grep "Appendix B" /tmp/sam-repo/stateless-software-engineering-framework.md | wc -l` — should be ≥4 (1 heading + ≥3 cross-references)
4. PR link confirmed

---

## Critical Constraints

- Do NOT invent SSF section numbers — use `[TODO: confirm section N.N]` until actual SSF structure is read
- Do NOT commit cross-repo changes without reading the existing SSF document style first — match heading levels, code block style, and table format
- The Appendix MUST be self-contained — a reader should not need to have read the rest of SSF to understand and apply the error recovery procedure
- All git commands in the procedure must be tested against a real scenario before inclusion — no theoretical commands
- T4 requires write access to `bitflight-devops/stateless-agent-methodology` — verify access before starting T4
- Conventional commit scope for external repo changes: `docs(ssf): ...`
- Back-reference: this plan file (`plan/tasks-5-sam-error-recovery.md`) must be mentioned in the PR body

---

## Next Steps (after plan completion)

```text
Backlog item "SAM: Error Recovery / Rollback Procedures" is now planned.

- Plan file: plan/tasks-5-sam-error-recovery.md
- To execute T1+T2: research and design tasks (no external repo needed)
- To execute T3: /python3-development:implement-feature sam-error-recovery (writes appendix draft)
- To execute T4: clone bitflight-devops/stateless-agent-methodology, apply draft
- To close when done: /work-backlog-item close SAM: Error Recovery / Rollback Procedures
```
