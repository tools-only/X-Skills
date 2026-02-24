# Architecture Spec: Fix CHANGELOG References to Nonexistent Files

## Document Metadata

- **Generated**: 2026-02-23
- **Feature**: conventional-commits-changelog-refs
- **Type**: Documentation remediation (no code changes)
- **Source**: `plan/feature-context-conventional-commits-changelog-refs.md`

---

## Executive Summary

This is a documentation-only fix. No Python code, CLI commands, or architectural components are being added. The goal is to remove or correct dead references in `plugins/conventional-commits/CHANGELOG.md` and `plugins/conventional-commits/skills/conventional-commits/SKILL.md`.

---

## Scope

### In Scope

1. **CHANGELOG.md** (`plugins/conventional-commits/CHANGELOG.md`)
   - Line 27: Reference to `docs/examples.md` (file does not exist)
   - Line 31: Placeholder release URL `https://github.com/owner/conventional-commits/releases/tag/v1.0.0`

2. **SKILL.md Related Skills** (`plugins/conventional-commits/skills/conventional-commits/SKILL.md`)
   - Lines 494-501: Related Skills section references commitlint, pre-commit, git-commit-helper, semantic-release
   - commitlint: EXISTS (`plugins/commitlint`)
   - pre-commit: EXISTS (`plugins/python3-development/skills/pre-commit`)
   - git-commit-helper: Skill is actually `commit-staged` (`.claude/skills/commit-staged`)
   - semantic-release: Does NOT exist as a skill

### Out of Scope

- Creating new skills (e.g., semantic-release)
- Modifying plugin.json
- Adding new documentation beyond fixing references

---

## Design Decisions (Resolved from Feature Context Questions)

### Decision 1: docs/examples.md

**Resolution**: Remove the CHANGELOG line claiming docs/examples.md exists.

**Rationale**: The file was never created. Creating it would require content authoring beyond the scope of "fix dead references." Removing the claim is honest and aligns with RT-ICA goal of removing dead references.

**Action**: Delete or rewrite line 27 to remove the "Added docs/examples.md" claim. Option: "Added extensive examples documentation covering common scenarios" (generic, no file path).

### Decision 2: Release URL

**Resolution**: Remove the placeholder link or replace with repository root releases URL.

**Rationale**: `owner/conventional-commits` is a placeholder. The conventional-commits plugin lives in `Jamie-BitFlight/claude_skills`. Plugin-specific releases may not exist. Per Keep a Changelog, the link is optional when no release exists.

**Action**: Remove the `[1.0.0]: https://github.com/owner/conventional-commits/releases/tag/v1.0.0` line, or replace with `https://github.com/Jamie-BitFlight/claude_skills/releases` if that is the canonical release location. Prefer removal if no v1.0.0 release exists for this plugin.

### Decision 3: Related Skills

**Resolution**: Update to reference only existing skills with correct names.

**Rationale**: Dead references cause activation failures. commit-staged is the actual skill name for "generate commit messages from diffs."

**Action**:
- commitlint: Keep (exists)
- pre-commit: Keep (exists; use activation syntax `Skill(command: "python3-development:pre-commit")` or equivalent)
- git-commit-helper: Replace with `commit-staged` — `Skill(command: "commit-staged")` or document as commit-staged
- semantic-release: Remove (no skill exists)

---

## File Change Specification

### CHANGELOG.md

| Location | Current | Target |
|----------|---------|--------|
| Line 27 | "Added docs/examples.md with 15 real-world usage scenarios" | "Added extensive examples documentation covering common scenarios" (or remove if redundant with line 21) |
| Line 31 | `[1.0.0]: https://github.com/owner/conventional-commits/releases/tag/v1.0.0` | Remove, or replace with valid URL |

### SKILL.md (Related Skills section)

| Reference | Current | Target |
|-----------|---------|--------|
| commitlint | `commitlint` | Keep; ensure activation syntax correct |
| pre-commit | `pre-commit` | Keep; ensure activation syntax correct |
| git-commit-helper | `git-commit-helper` | Change to `commit-staged` |
| semantic-release | `semantic-release` | Remove |

---

## Verification Requirements

1. **CHANGELOG**: No references to non-existent files; no placeholder URLs
2. **SKILL.md**: All Related Skills exist and use correct activation syntax
3. **Plugin validation**: `claude plugin validate plugins/conventional-commits/` passes
4. **plugin_validator**: Run if available; no new errors

---

## Success Criteria

- [ ] CHANGELOG.md contains no dead file references
- [ ] CHANGELOG.md contains no placeholder release URLs (or has valid URL)
- [ ] SKILL.md Related Skills section references only existing skills
- [ ] Plugin validation passes
