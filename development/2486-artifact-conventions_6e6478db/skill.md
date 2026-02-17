# Artifact Conventions

SAM artifact naming, file layout, and cross-referencing conventions for the development harness.

---

## Principle

Every stage in the SAM pipeline produces a file-based artifact. These artifacts are the only communication channel between stages. No stage relies on conversation memory — each reads its predecessor's artifact and writes its own. This is the stateless property of SAM.

---

## Cross-Reference Token Pattern

**Format:** `ARTIFACT:{TYPE}({SCOPE_OR_ID})`

**Types:**

- `ARTIFACT:DISCOVERY({feature-slug})` — S1 output
- `ARTIFACT:PLAN({feature-slug})` — S2 output
- `ARTIFACT:CONTEXT({feature-slug})` — S3 output (amends S2)
- `ARTIFACT:TASK({task-id})` — S4 output (one per task)
- `ARTIFACT:EXECUTION({task-id})` — S5 output (one per task)
- `ARTIFACT:REVIEW({feature-slug})` — S6 output
- `ARTIFACT:VERIFICATION({feature-slug})` — S7 output

**Usage in artifacts:** Each artifact includes a header block with references to predecessor and successor artifacts.

```markdown
---
artifact: ARTIFACT:PLAN(add-jwt-auth)
predecessor: ARTIFACT:DISCOVERY(add-jwt-auth)
successor: ARTIFACT:CONTEXT(add-jwt-auth)
feature: add-jwt-auth
stage: S2
created: 2026-02-15
---
```

---

## File Layout

All artifacts are stored in `.planning/harness/` relative to the project root.

**Directory structure:**

```text
.planning/
  harness/
    discovery-{feature-slug}.md
    plan-{feature-slug}.md
    context-{feature-slug}.md
    task-{task-id}-{task-slug}.md     (one per task)
    execution-{task-id}-{task-slug}.md (one per task)
    review-{feature-slug}.md
    verification-{feature-slug}.md
```

---

## File Naming Conventions

### Feature-Level Artifacts

**Pattern:** `{stage-prefix}-{feature-slug}.md`

**Stage prefixes by stage:**

- S1 — `discovery`
- S2 — `plan`
- S3 — `context`
- S6 — `review`
- S7 — `verification`

**Feature slug rules:**

- Derived from the feature name or request
- Lowercase, hyphen-separated
- Max 50 characters
- No special characters beyond hyphens

**Examples:**

- `discovery-add-jwt-auth.md`
- `plan-refactor-payment-module.md`
- `review-add-csv-export.md`

### Task-Level Artifacts

**Pattern:** `{stage-prefix}-{task-id}-{task-slug}.md`

**Stage prefixes for task-level:**

- S4 — `task`
- S5 — `execution`

**Task ID rules:**

- Three-digit zero-padded integer
- Assigned sequentially during S4 Task Decomposition
- IDs are stable — not renumbered after creation

**Task slug rules:**

- Brief description of the task
- Same formatting rules as feature slug

**Examples:**

- `task-001-add-jwt-middleware.md`
- `task-002-add-token-validation.md`
- `execution-001-add-jwt-middleware.md`
- `execution-002-add-token-validation.md`

---

## Artifact Types

### DISCOVERY

**Stage:** S1
**Purpose:** Capture requirements, codebase context, constraints, and resolved role assignments.

**Required sections:**

- Feature request (original, verbatim)
- Codebase context (structure, relevant files, patterns)
- Constraints identified (bound and unbound)
- Role resolution results (which agents assigned to which roles)
- Quality gates configured
- ARL touchpoint assessment for Gate 1

### PLAN

**Stage:** S2
**Purpose:** Development plan with task graph, acceptance criteria, and RT-ICA analysis.

**Required sections:**

- RT-ICA gap analysis (present, partial, missing information)
- Task graph with dependencies
- Task skeletons with acceptance criteria
- Quality gate schedule (which gates run at which stages)
- Unblock paths for tasks with missing information

### CONTEXT

**Stage:** S3
**Purpose:** Plan validated against actual codebase state.

**Required sections:**

- Validation results per plan element
- Codebase discrepancies found (if any)
- Plan amendments (if needed)
- Verified integration points
- Confirmed dependency availability

### TASK

**Stage:** S4
**Purpose:** Individual executable task file.

**Required sections:**

- Task ID and description
- Predecessor artifact reference
- Assigned agent (from role resolution)
- Inputs (files to read, context needed)
- Acceptance criteria (testable, specific)
- Quality gates to run after completion
- Dependencies on other tasks (if any)

### EXECUTION

**Stage:** S5
**Purpose:** Record of what was implemented for a task.

**Required sections:**

- Task reference
- Files created or modified (with paths)
- Implementation decisions and rationale
- Quality gate results (pass/fail per gate)
- Deviations from plan (if any)

### REVIEW

**Stage:** S6
**Purpose:** Forensic review comparing execution against acceptance criteria.

**Required sections:**

- Per-task verdict (COMPLETE or NEEDS_WORK)
- Acceptance criteria evaluation (pass/fail per criterion)
- Quality gate summary
- Regression check results
- NEEDS_WORK details (what failed, what to fix)

### VERIFICATION

**Stage:** S7
**Purpose:** Final certification against original requirements.

**Required sections:**

- Per-requirement verdict (MET or NOT_MET)
- Traceability matrix (requirement to task to implementation)
- Full quality gate results
- Integration verification results
- Final verdict (CERTIFIED or NOT_CERTIFIED)
- NOT_CERTIFIED details (what requirements are unmet)

---

## Cross-Referencing Between Artifacts

Each artifact references its immediate predecessor and the tasks it relates to.

**Forward references:** When S4 creates tasks, each task artifact includes `predecessor: ARTIFACT:PLAN({feature-slug})`.

**Backward references:** When S6 reviews tasks, the review artifact includes references to all `ARTIFACT:EXECUTION({task-id})` artifacts.

**Traceability chain:** S7 verification builds a complete traceability matrix linking requirements (from S1) through plan (S2), tasks (S4), execution (S5), and review (S6).

---

## Coexistence with Other Planning Tools

The `.planning/harness/` directory is scoped to the development harness. Other tools use their own subdirectories:

- `.planning/gsd/` — GSD (Get Stuff Done) planning tool
- `.planning/backlog/` — Backlog items and grooming
- `.planning/` root — Shared or tool-agnostic planning documents

The harness never reads or writes outside `.planning/harness/`. No naming collisions occur because each tool uses its own prefix and subdirectory.

---

## Cleanup

Artifacts persist after feature completion for auditability. The harness does not auto-delete artifacts.

**Manual cleanup:** Users can delete `.planning/harness/` or specific feature artifacts after verification.

**Gitignore recommendation:** Add `.planning/` to `.gitignore` unless the team wants planning artifacts in version control.

---

## Sources

- SAM methodology: <https://github.com/bitflight-devops/stateless-agent-methodology>
- Default development flow: [./default-development-flow.md](./default-development-flow.md)
