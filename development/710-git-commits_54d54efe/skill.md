---
trigger: always_on
glob:
description: Git commit rules and quality gate protections
---

# Git Commit Rules

## Quality Gate Protection

**NEVER use `--no-verify`, `--skip-hooks`, or any flag that bypasses git hooks or pre-commit checks.**

These quality gates exist for a reason. When a hook fails:

1. Read the error message carefully
2. Check the relevant configuration files (e.g., `.pre-commit-config.yaml`)
3. Consult reference documentation in this repository (e.g., skills, plugins)
4. Fix the underlying issue properly

Only the user may choose to bypass these protections—never the AI assistant.

## Conventional Commits

This repository enforces [Conventional Commits](https://www.conventionalcommits.org/) via `prek` (a Rust-based drop-in replacement for pre-commit).

### Configuration

The prek hook uses `--strict --force-scope`, meaning:

- A **scope is required** in all commit messages
- Format: `<type>(<scope>): <description>`

### Valid Commit Types

| Type       | Description                         |
| ---------- | ----------------------------------- |
| `feat`     | New feature                         |
| `fix`      | Bug fix                             |
| `docs`     | Documentation only                  |
| `style`    | Code style changes (formatting)     |
| `refactor` | Code change with no bug fix/feature |
| `perf`     | Performance improvement             |
| `test`     | Adding or correcting tests          |
| `build`    | Build system or dependency changes  |
| `ci`       | CI configuration changes            |
| `chore`    | Maintenance (no src/test changes)   |
| `revert`   | Reverts a previous commit           |

### Examples

```text
feat(auth): add user authentication
fix(parser): handle null pointer exception
docs(claude): update skills and add refactor plan
refactor(api): extract validation logic
ci(pipeline): add Node 18 to test matrix
```

### Resources

- Skill reference: `plugins/conventional-commits/skills/conventional-commits/SKILL.md`
- Pre-commit config: `.pre-commit-config.yaml`
