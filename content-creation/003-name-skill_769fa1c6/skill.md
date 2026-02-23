---
name: update-readme
description: >
  This skill should be used when the user asks to "update my README",
  "refresh the README", "README is outdated", "sync README with the codebase",
  "improve my README", "keep README up to date", "fix the README",
  "my README is stale", "update the readme file", "README needs updating",
  "update readme from git history", or "readme is out of date".
allowed-tools:
  - Read
  - Glob
  - Grep
  - Write
  - Edit
  - Task
  - AskUserQuestion
  - Skill
---

# Update README

Refresh an existing `README.md` using current codebase state, git history, and
changelog intelligence. Always updates `CHANGELOG.md` first so changelog content
is available when revising README sections.

## Step 1 — Check README Existence

Glob for `README.md`, `readme.md`, and `Readme.md` in the working directory.

If not found, or if the file has fewer than 5 lines of real content:
- Inform the user no substantial README was found
- Invoke `Skill: claude-coding:make-readme`
- Stop — make-readme handles creation from scratch

If a substantial README exists: continue.

## Step 2 — Update Changelog

Invoke `Skill: claude-coding:make-changelog`.

make-changelog detects fill mode automatically — it adds only versions missing
since the last documented entry. If the repo has no git history, make-changelog
will report an error; note it and continue without changelog context.

## Step 3 — Parallel Research

Record the current README line count before edits.

Launch the three agents below **in a single message** (all in parallel):

**Agent A — README Audit** (`subagent_type: Explore`)

Prompt:
```
Read README.md. For each section (every H1/H2/H3 heading), record:
- Section name and approximate line count
- Stale content: version numbers, paths, commands, or features that look incorrect
- Placeholder text never replaced (e.g. "YOUR_USERNAME", "TODO", angle-bracket values)
- Thin sections with fewer than 3 lines of real content

Also identify which standard sections are MISSING:
installation, usage, configuration, API reference, contributing, license.

Return structured notes:
- existing_sections: list of heading names
- stale_items: list of {location, description} objects
- missing_sections: list
- thin_sections: list
```

**Agent B — Codebase Scan** (`subagent_type: Explore`)

Prompt:
```
Scan the working directory to extract:
- Project name and current version from the first found manifest:
  package.json, pyproject.toml, Cargo.toml, go.mod, .claude-plugin/plugin.json
- Primary language and framework (file extensions, config files)
- Top-level directory structure, max 2 levels deep
- Key public-facing files: entry points, config templates, example files, CLI scripts
- Documentation files beyond README: CONTRIBUTING.md, docs/, ARCHITECTURE.md, etc.

Return: project_name, version, language, framework, structure_summary,
key_files (paths), extra_docs (paths).
```

**Agent C — Git History Since Last README Touch** (`subagent_type: Explore`)

Prompt:
```
Run: git log --follow -1 --format="%ai" -- README.md
This gives the date README.md was last committed.
If README.md has no git history, use the initial commit date:
  git log --reverse --format="%ai" | head -1

Then run: git log --since="[date from above]" --format="%s" --no-merges

Categorize each commit subject into:
- new_features: new capabilities users can invoke
- breaking_changes: removed or incompatible changes
- deprecations: features flagged for future removal
- significant_fixes: user-visible bug fixes
- other: everything else

Skip: merge commits, CI/CD changes, formatting/whitespace, version-bump commits.

Return: readme_last_updated (ISO date), changes_since (object with 5 category lists).
```

Wait for all three agents to complete before Step 4.

## Step 4 — Read Updated Changelog

Read `CHANGELOG.md` and extract the most recent version section (the first `## [x.y.z]`
block). This provides structured change context for README sections like Features.

## Step 5 — Synthesize and Update

Apply targeted updates to `README.md`. Preserve existing style, header format,
badge patterns, and overall structure — do not restructure unless clearly broken.

Apply in this priority order — highest-confidence fixes first, judgment-heavy
additions last, because objective errors (wrong version, broken path) have
deterministic answers while structural additions require codebase context to get right:

1. **Version/badge updates** — update all version numbers and shields.io badge
   URLs to the version from Agent B
2. **Stale content** — correct outdated paths, commands, flags, or removed features
   flagged by Agent A
3. **Placeholder cleanup** — replace any unfilled placeholders with real values
   from Agent B
4. **New features** — add or expand the Features section using Agent C's
   `new_features` list and the latest CHANGELOG section
5. **Missing sections** — add any standard sections identified by Agent A as
   absent, with minimal viable content drawn from Agent B's codebase scan
6. **Thin sections** — expand underdeveloped sections using Agent B context

If no updates are needed after applying all criteria, skip edits and proceed to Step 6.

For each change, make the smallest edit that achieves accuracy. Use `Edit` for
targeted section updates — it preserves git blame and is safer for incremental
changes. Use `Write` only if more than 60% of the file changes, because at that
point the document is being regenerated, not updated.

## Step 6 — Report

Output a diff summary after writing:
- Line count: N → M lines
- Sections updated: list with one-line description of each change
- Sections added: list
- Sections removed: list with reason (if any)
- Version bumped: old → new (if applicable)
- Changelog: versions added by make-changelog
