---
name: updateclaudemd
description: >
  This skill should be used when the user asks to update, refresh, optimize, or
  clean up CLAUDE.md, or when documentation is stale, verbose, or out of sync
  with codebase reality. Triggers on "update CLAUDE.md", "refresh the docs",
  "sync CLAUDE.md with the codebase", "optimize project instructions",
  "clean up CLAUDE.md", "improve CLAUDE.md", "fix CLAUDE.md".
allowed-tools:
  - Read
  - Glob
  - Grep
  - Edit
  - Write
  - Bash(wc:*)
  - AskUserQuestion
hooks:
  PreToolUse:
    - matcher: "Write"
      command: "cp CLAUDE.md CLAUDE.md.bak 2>/dev/null || true"
      once: true
---

# Update and Optimize CLAUDE.md

## Principles

Content in CLAUDE.md is justified only if it changes how Claude acts in the next session. Apply this test to every line — remove anything that does not influence behavior. Aim for 150-250 lines; favor patterns and principles over verbose documentation.

## Step 1: Read Current State and Explore Codebase

If no `CLAUDE.md` exists, inform the user and begin from scratch using codebase exploration only. Otherwise, read the existing `CLAUDE.md` and note its line count, section structure, and any content that looks stale or duplicated. Then explore the project systematically: read configuration files (package.json, Cargo.toml, requirements.txt, etc.), map the directory structure, identify the tech stack, and note established patterns and conventions. Focus on what would change how a future Claude session works in this codebase.

Stop exploration when the major directories, tech stack, and 3-5 key patterns are identified — patterns that affect how code should be written, tested, or deployed in this project. Do not attempt exhaustive coverage.

## Step 2: Reconcile and Optimize

Compare codebase reality with current documentation. Apply the governing principle to every line.

Verify that commands, paths, dependencies, and environment variables are still accurate.

The dividing line: remove content Claude would infer from reading the codebase itself; keep content that would cause errors or suboptimal behavior without explicit documentation.

Cut: code blocks duplicating source files, verbose explanations, one-time setup troubleshooting, philosophical "why this matters" sections, anything Claude can derive from file inspection.

Keep: architecture decisions, essential development commands, coding conventions, critical gotchas, non-obvious behaviors that need documentation to prevent errors.

Proceed to Step 3 when every section has been evaluated against the governing principle and all inaccurate references are corrected.

## Step 3: Structure for Scanning

Organize into scannable sections appropriate for this project type. Common sections: project overview, development commands, architecture principles, project structure, coding conventions, and development notes. Add conditional sections (design system, database schema, environment variables) only when the project warrants them.

If the existing CLAUDE.md has content whose purpose is unclear, use AskUserQuestion to confirm before removing it.

## Step 4: Write the Optimized Version

Write the updated CLAUDE.md. Ensure scannable headers, concise actionable content, accurate reflection of the current codebase, and minimal redundancy.

A `CLAUDE.md.bak` backup was created automatically before this write. Inform the user it exists for diffing and suggest removing it after review.

## Output Summary

Summarize what changed: line count before vs after, sections added or removed, and key corrections applied.
