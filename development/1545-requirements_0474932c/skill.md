# Requirements: Pre-Execution Validation for Plan and Design Commands

## Metadata
- **Feature**: plan-validation-checks
- **Status**: APPROVED
- **Created**: 2026-02-05
- **Author**: Claude Opus 4.5

---

## 1. Problem Statement

ZERG commands `/z:plan` and `/z:design` can execute redundantly when:
- The objective has already been implemented in recent commits
- An open PR already addresses the feature
- Codebase already contains matching implementations

This wastes time and can create conflicts with in-flight work.

---

## 2. Functional Requirements

### FR-1: Pre-Execution Validation Sequence

Before executing either `/z:plan` or `/z:design`, perform this 5-step validation:

1. **Extract Objective** — Read the plan/requirements file and extract stated objective
2. **Check Recent Commits** — Run `git log --oneline -20` to review recent work
3. **Check Open PRs** — Run `gh pr list` to identify in-flight work
4. **Search Codebase** — Grep for implementations matching the plan's targets
5. **Conflict Detection** — If objective appears completed or conflicts exist, STOP

### FR-2: Conflict Resolution Options

When validation fails, present:
- `Update plan` — Modify the spec to reflect current state
- `Archive it` — Move spec to `.gsd/specs/_archived/`
- `Proceed anyway` — Continue with explicit user override

### FR-3: Validation Pass Criteria

Proceed only if:
- No commits in last 20 match the feature name or key terms
- No open PRs match the feature name or key terms
- Grep for key identifiers returns <5 matches (configurable threshold)

---

## 3. Affected Files

| File | Change |
|------|--------|
| `zerg/data/commands/plan.core.md` | Add Phase 0: Validation before Phase 1 |
| `zerg/data/commands/design.core.md` | Add Phase 0: Validation before Load Context |

---

## 4. Implementation Details

### 4.1 New Section in plan.core.md (Insert before "## Enter Plan Mode")

```markdown
## Phase 0: Pre-Execution Validation

Before proceeding, validate this plan hasn't been superseded:

1. **Extract Objective**
   - Read `.gsd/specs/$FEATURE/requirements.md` if exists
   - Identify key terms: feature name, main components, file patterns

2. **Check Recent Commits**
   ```bash
   git log --oneline -20 | grep -i "$FEATURE"
   ```
   - Flag if any commits mention the feature name

3. **Check Open PRs**
   ```bash
   gh pr list --state open | grep -i "$FEATURE"
   ```
   - Flag if any open PRs match

4. **Search Codebase**
   - Grep for key implementation patterns from the requirements
   - Flag if substantial matches found (>5 files)

5. **Validation Decision**
   IF any checks flag potential conflicts:
     STOP and present:
     ```
     ⚠️  VALIDATION WARNING

     Potential conflict detected:
     - [Commits/PRs/Code] matching "{feature}" found

     Options:
     1. Update plan - Revise spec to account for existing work
     2. Archive - Move to .gsd/specs/_archived/
     3. Proceed anyway - Override and continue
     ```

     Use AskUserQuestion to get user decision.

   IF validation passes:
     Continue to Phase 1.
```

### 4.2 New Section in design.core.md (Insert before "## Load Context")

```markdown
## Phase 0: Pre-Execution Validation

Before proceeding, validate this design is still needed:

1. **Read Requirements**
   - Load `.gsd/specs/$FEATURE/requirements.md`
   - Extract key objectives and target files

2. **Check Recent Commits**
   ```bash
   git log --oneline -20
   ```
   - Compare commit messages against requirements objectives

3. **Check Open PRs**
   ```bash
   gh pr list --state open
   ```
   - Check for PRs implementing similar features

4. **Grep Targets**
   - For each file in requirements' "Files to Create/Modify":
     ```bash
     ls {target_file} 2>/dev/null && echo "EXISTS: {target_file}"
     ```
   - For key function/class names mentioned:
     ```bash
     grep -r "{identifier}" zerg/ tests/ --include="*.py" -l | head -5
     ```

5. **Validation Decision**
   IF requirements.md missing:
     ERROR: Run /z:plan first

   IF target files already exist OR key identifiers found:
     STOP and present conflict resolution options

   IF validation passes:
     Continue to Load Context
```

---

## 5. Non-Functional Requirements

- **NFR-1**: Validation adds <5 seconds to command startup
- **NFR-2**: Validation is skippable with `--skip-validation` flag
- **NFR-3**: Archived specs retain full history in git

---

## 6. Acceptance Criteria

- [ ] `/z:plan feature-x` checks git log before starting
- [ ] `/z:plan feature-x` checks gh pr list before starting
- [ ] `/z:design` checks if target files exist before designing
- [ ] Conflict detection presents clear options to user
- [ ] `--skip-validation` bypasses all checks
- [ ] Archived specs go to `.gsd/specs/_archived/{feature}/`

---

## 7. Out of Scope

- Automatic conflict resolution (requires user decision)
- Semantic analysis of code changes (just string matching)
- Cross-repository validation

---

## 8. Open Questions

None — requirements are self-contained from user request.
