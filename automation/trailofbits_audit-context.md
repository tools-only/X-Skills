---
name: trailofbits:audit-context
source: https://raw.githubusercontent.com/trailofbits/skills/main/plugins/audit-context-building/commands/audit-context.md
original_path: plugins/audit-context-building/commands/audit-context.md
source_repo: trailofbits/skills
category: automation
subcategory: workflow
tags: ['automation']
collected_at: 2026-02-01T04:15:15.908492
file_hash: 16cc3a488df3eae72b77e3ad55102acabc09320fd1e65c5520fd72ec6ece4b33
---

---
name: trailofbits:audit-context
description: Builds deep architectural context before vulnerability hunting
argument-hint: "<codebase-path> [--focus <module>]"
allowed-tools:
  - Read
  - Grep
  - Glob
  - Bash
  - Task
---

# Build Audit Context

**Arguments:** $ARGUMENTS

Parse arguments:
1. **Codebase path** (required): Path to codebase to analyze
2. **Focus** (optional): `--focus <module>` for specific module analysis

Invoke the `audit-context-building` skill with these arguments for the full workflow.
