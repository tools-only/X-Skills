---
paths:
  - "**/SKILL.md"
  - "**/agents/**/*.md"
  - "**/commands/**/*.md"
---

# Frontmatter Requirements

## Skills

- `name`: Required — lowercase, hyphens, must match directory name, satisfies `^[a-z][a-z0-9-]*$`
- `description`: Optional (uses first paragraph if omitted)
- `tools`: Must be comma-separated string — `Read, Grep, Glob` — not a YAML array

## Agents

- `name`: Required — lowercase, hyphens, max 64 chars
- `description`: Required — include trigger keywords, max 1024 chars, no colons except in URLs
- `model`: Must be `sonnet`, `opus`, `haiku`, or `inherit` if specified
- `tools`: Must be comma-separated string (not YAML array)
- No YAML multiline indicators (`>-`, `|-`, `>`, `|`) in any field

## Commands

- `description`: Required
- `allowed-tools`: Must be comma-separated string (not YAML array)

## Validator Auto-Fix

Run after writing or editing any frontmatter file:

```bash
uv run plugins/plugin-creator/scripts/plugin_validator.py --fix {path}
```

The validator auto-adds `name:` derived from the directory name when absent (plugin skills only).

**SOURCE:** `plugin_validator.py` and [agentskills.io specification](https://agentskills.io/specification)
