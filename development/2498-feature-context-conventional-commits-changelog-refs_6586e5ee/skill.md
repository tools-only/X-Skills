# Feature Context: Fix CHANGELOG References to Nonexistent Files

## Document Metadata

- **Generated**: 2026-02-23
- **Input Type**: simple_description (backlog item)
- **Source**: Plugin code review session 2026-02-21, BACKLOG.md
- **Status**: DISCOVERY_COMPLETE

---

## Original Request

CHANGELOG references files that do not exist in the repository. Additionally, related skills referenced in the plugin do not exist. All dead references need to be either created or removed.

**Files**: `plugins/conventional-commits/` (CHANGELOG and skill cross-references)

---

## Core Intent Analysis

### WHO (Target Users)

- Plugin maintainers and contributors
- Users reading CHANGELOG for release history
- AI agents following skill cross-references

### WHAT (Desired Outcome)

- CHANGELOG contains only references to files that exist
- Release URL in CHANGELOG points to a valid location (or is removed if no releases exist)
- Related Skills section in SKILL.md references only skills that exist in the repository

### WHEN (Trigger Conditions)

- Plugin validation runs (`claude plugin validate`)
- User or agent follows CHANGELOG links
- User or agent activates a Related Skill from conventional-commits SKILL.md

### WHY (Problem Being Solved)

- Dead references cause confusion and failed lookups
- Placeholder URLs (e.g., `github.com/owner/conventional-commits`) mislead users
- Non-existent skill references waste agent context and produce errors

---

## Codebase Research

### Similar Patterns Found

#### Pattern 1: CHANGELOG with valid release links

- **Location**: `plugins/clang-format/CHANGELOG.md`
- **Relevance**: Other plugins in the repo have CHANGELOGs; can compare structure and link patterns
- **Reusable**: Keep a Changelog format, fix link targets

#### Pattern 2: Plugin skill cross-references

- **Location**: `plugins/commitlint/skills/commitlint/SKILL.md` — references conventional-commits
- **Relevance**: Skills reference each other; need consistent naming (commit-staged vs git-commit-helper)
- **Reusable**: Use canonical skill names from `.claude/skills/` and `plugins/*/skills/`

#### Pattern 3: docs/ folder usage in plugins

- **Location**: `plugins/conventional-commits/` — no docs/ folder exists
- **Relevance**: CHANGELOG line 27 claims "Added docs/examples.md" but docs/ does not exist
- **Reusable**: Either create docs/examples.md or remove the CHANGELOG claim

### Existing Infrastructure

- `plugins/conventional-commits/CHANGELOG.md` — lines 27-31 contain dead references
- `plugins/conventional-commits/skills/conventional-commits/SKILL.md` — Related Skills section (lines 494-501)
- `plugins/commitlint/` — exists; skill at `plugins/commitlint/skills/commitlint/SKILL.md`
- `plugins/python3-development/skills/pre-commit/` — exists
- `.claude/skills/commit-staged/SKILL.md` — exists (generates conventional commit messages; referenced elsewhere as "git-commit-helper")
- No `plugins/semantic-release/` or `.claude/skills/semantic-release/` — semantic-release skill does not exist

### Code References

- `plugins/conventional-commits/CHANGELOG.md:27` — "Added docs/examples.md with 15 real-world usage scenarios"
- `plugins/conventional-commits/CHANGELOG.md:31` — `[1.0.0]: https://github.com/owner/conventional-commits/releases/tag/v1.0.0` (placeholder)
- `plugins/conventional-commits/skills/conventional-commits/SKILL.md:497-500` — Related Skills: commitlint, pre-commit, git-commit-helper, semantic-release

---

## Use Scenarios

### Scenario 1: User follows CHANGELOG link to docs/examples.md

**Actor**: Developer reading release notes
**Trigger**: Clicks or navigates to docs/examples.md referenced in CHANGELOG
**Goal**: View 15 real-world usage scenarios
**Expected Outcome**: File exists and displays content
**Current State**: 404 / file not found

### Scenario 2: User follows release URL

**Actor**: Developer checking release history
**Trigger**: Clicks `[1.0.0]` link in CHANGELOG
**Goal**: View GitHub release for v1.0.0
**Expected Outcome**: Valid release page or no misleading placeholder
**Current State**: Placeholder URL points to non-existent repo (owner/conventional-commits)

### Scenario 3: Agent activates Related Skill

**Actor**: AI agent using conventional-commits skill
**Trigger**: Skill suggests "git-commit-helper" or "semantic-release" for related tasks
**Goal**: Activate the referenced skill
**Expected Outcome**: Skill loads successfully
**Current State**: git-commit-helper may map to commit-staged; semantic-release skill does not exist

---

## Gap Analysis

### Identified Gaps

| # | Category | Gap Description | Impact |
|---|----------|-----------------|--------|
| 1 | Scope | docs/examples.md: Create vs remove? | If create: need content. If remove: update CHANGELOG text |
| 2 | Scope | Release URL: Fix to real repo vs remove? | Repo is Jamie-BitFlight/claude_skills; conventional-commits is a plugin, not a separate repo |
| 3 | Integration | git-commit-helper: Alias for commit-staged? | .claude/skills/README uses "git-commit-helper" but skill is commit-staged |
| 4 | Integration | semantic-release: Remove reference or create skill? | No semantic-release skill exists |

---

## Questions Requiring Resolution

### Q1: docs/examples.md — create or remove reference?

- **Category**: Scope
- **Gap**: CHANGELOG claims "Added docs/examples.md with 15 real-world usage scenarios" but file does not exist
- **Question**: Should we (A) create docs/examples.md with content, or (B) remove the CHANGELOG line claiming it exists?
- **Options**:
  - A) Create `plugins/conventional-commits/docs/examples.md` with 15 usage scenarios (aligns with CHANGELOG claim)
  - B) Remove the line from CHANGELOG (honest: we never created it)
- **Why It Matters**: Creating requires content authoring; removing is faster but changes historical record
- **Resolution**: _pending_

### Q2: Release URL — fix or remove?

- **Category**: Scope
- **Gap**: `[1.0.0]: https://github.com/owner/conventional-commits/releases/tag/v1.0.0` is placeholder
- **Question**: Should we (A) point to Jamie-BitFlight/claude_skills releases (if any), (B) remove the link, or (C) use a different format?
- **Options**:
  - A) Use `https://github.com/Jamie-BitFlight/claude_skills/releases` if releases exist for this plugin
  - B) Remove the link line (Keep a Changelog allows this when no release URL exists)
  - C) Use `./` or relative link if applicable
- **Why It Matters**: Placeholder misleads; real URL improves navigation
- **Resolution**: _pending_

### Q3: Related Skills — fix names and remove dead refs?

- **Category**: Integration
- **Gap**: git-commit-helper and semantic-release may not exist as skills
- **Question**: Update Related Skills to reference only existing skills?
- **Options**:
  - A) commitlint → exists (`plugins/commitlint`). pre-commit → exists (`plugins/python3-development/skills/pre-commit`). git-commit-helper → change to commit-staged (`.claude/skills/commit-staged`). semantic-release → remove (no skill)
  - B) Keep git-commit-helper if it's a documented alias; remove semantic-release
- **Why It Matters**: Dead references cause activation failures
- **Resolution**: _pending_

---

## Goals (Pending Resolution)

_These goals will be finalized after questions are resolved._

1. Remove or fix all dead references in `plugins/conventional-commits/CHANGELOG.md`
2. Remove or fix dead references in `plugins/conventional-commits/skills/conventional-commits/SKILL.md` Related Skills section
3. Run `claude plugin validate plugins/conventional-commits/` and ensure no reference-related failures
4. Optionally run `plugin_validator` if applicable

---

## Next Steps

After questions are resolved:

1. Update "Resolution" fields in Questions section
2. Finalize Goals section
3. Proceed to architecture design (minimal for doc-only fix)
4. Create task decomposition
5. Execute implementation
