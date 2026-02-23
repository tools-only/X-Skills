---
name: make-changelog
description: >
  This skill should be used when the user asks to "create a changelog",
  "generate a changelog", "update my changelog", "fill in the changelog",
  "add a changelog", "CHANGELOG is missing entries", "changelog is out of date",
  "what's missing from my changelog", "changelog from git history", "write
  changelog", "release notes", or says "my project needs a CHANGELOG".
model: sonnet
context: fork
agent: general-purpose
allowed-tools:
  - Bash(git:*)
  - Bash(python3:*)
  - Read
  - Write
  - Edit
  - Task
  - AskUserQuestion
---

# Changelog Generator

Create or update `CHANGELOG.md` from git history using Keep-a-Changelog format.
Launch one haiku subagent per version range for parallel, token-efficient processing.

## Step 1 — Assess Project State

Run the range-planning script to validate the repo and gather all version ranges:

```bash
python3 ${CLAUDE_PLUGIN_ROOT}/skills/make-changelog/scripts/list_ranges.py --output json
```

- Exit 2 = not a git repository, or repository has no commits. Stop and report the
  error message from stderr to the user.
- Exit 0 = proceed. JSON output contains all ranges the skill would cover.

Also check for an existing changelog:

```bash
ls CHANGELOG.md CHANGELOG 2>/dev/null
```

If `CHANGELOG.md` exists, read it to identify the last documented version (the most
recent `## [x.y.z]` heading).

## Step 2 — Determine Scope

**No tags exist:** Script returns a single Unreleased range. Proceed in fresh mode.

**Tags exist, no CHANGELOG.md:** Fresh mode — process all ranges the script returned.

**Tags exist, CHANGELOG.md exists:** Fill mode — re-run the script with `--since-tag`
to get only the uncovered ranges:

```bash
python3 ${CLAUDE_PLUGIN_ROOT}/skills/make-changelog/scripts/list_ranges.py \
  --since-tag <last_documented_tag> --output json
```

If the intent is ambiguous (user said "update changelog" but the file has extensive
existing content), ask the user to choose:

- "Fill missing only" — add new versions/unreleased since last entry
- "Rebuild from scratch" — regenerate the full file from git history
- "Unreleased only" — refresh only the `[Unreleased]` section

## Step 3 — Plan Subagents

Inspect the JSON output from Step 1 or 2. Each range object has:
- `label` — version string (`1.2.0` or `Unreleased`)
- `from_ref` — lower git ref (empty string = from initial commit)
- `to_ref` — upper git ref
- `date` — tag date or today
- `commit_count` — number of commits in range

**Skip ranges where `commit_count == 0`** — no subagent needed.

**Cap at 12 subagents.** If the non-empty range list exceeds 12, ask the user whether
to limit to the most recent 12 or process all (noting it will take longer).

## Step 4 — Launch Haiku Subagents

Spawn one Task per non-empty range in a **single message** (parallel). Use
`model: "haiku"` for all subagents. Each subagent receives this self-contained
prompt (substitute values from the range object):

```
You are generating one section of a CHANGELOG.md.

Version: [label]
Date: [date]
Git command: git log [from_ref]..[to_ref] --format="%s"
(If from_ref is empty, use: git log [to_ref] --format="%s")

Categorize by user-observable impact, not by commit prefix. The conventional
commit prefix (feat:, fix:, etc.) is a hint only — a commit labelled "feat:"
that clearly corrects a bug belongs in Fixed, not Added. Use the prefix as a
starting signal and override it when the commit subject contradicts it.

Sections and classification signals:
- Added: new capabilities users can invoke (add, introduce, implement, support)
- Changed: modified behavior of existing features (update, change, refactor, improve)
- Deprecated: features explicitly flagged for future removal
- Removed: capabilities or endpoints deleted (remove, drop, delete)
- Fixed: bugs corrected, regardless of prefix (fix, resolve, patch, correct, hotfix)
- Security: vulnerability patches (security, CVE, sanitize, escape)

Skip: merge commits, CI/CD configuration changes, version-bump commits,
formatting-only changes.

Use present tense, imperative mood: "Add X", not "Added X" or "adds X".
Omit empty sections. Omit internal variable names and implementation details.

Output ONLY the markdown block, no preamble:

## [VERSION] - DATE

### Added
- ...

### Fixed
- ...
```

Collect all Task results before Step 5.

## Step 5 — Assemble Changelog

Order version blocks newest-first. Full document structure:

```markdown
# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

[content from unreleased subagent, or omit section if empty]

## [1.2.0] - 2024-03-15

[content from subagent]
```

**Fill mode:** Insert new version sections immediately after the `# Changelog` header
block, before the first existing `## [` entry. Preserve all existing content exactly.

**Fresh/rebuild mode:** Write the complete file.

## Step 6 — Write and Report

**Fresh mode:** Write `CHANGELOG.md`. Confirm overwrite with `AskUserQuestion` if
the file already exists.

**Fill mode:** Use `Edit` to insert new sections at the correct position.

After writing, report:
- Mode used (fresh / fill / unreleased-only)
- Versions covered (label, date, commit count — from script output)
- Number of subagents used
- Any ranges that were empty (skipped)
- Any commits that could not be confidently categorized (flag for user review)

## Scripts

`scripts/list_ranges.py` — queries git tags and computes version ranges. Invoked in
Steps 1 and 2. Exit 2 = repo invalid or empty; exit 0 = success.

```bash
# All ranges (text preview)
python3 ${CLAUDE_PLUGIN_ROOT}/skills/make-changelog/scripts/list_ranges.py

# All ranges (JSON for skill use)
python3 ${CLAUDE_PLUGIN_ROOT}/skills/make-changelog/scripts/list_ranges.py --output json

# Fill mode: only ranges newer than a specific tag
python3 ${CLAUDE_PLUGIN_ROOT}/skills/make-changelog/scripts/list_ranges.py \
  --since-tag v1.2.0 --output json
```
