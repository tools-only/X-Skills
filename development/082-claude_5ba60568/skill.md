# Plugin Creator Scripts - AI-Facing Documentation

## Purpose

This directory contains executable scripts for plugin development, validation, and maintenance.

## Scripts Overview

### auto_sync_manifests.py

Automatically maintains `plugin.json` and `marketplace.json` during pre-commit. Detects CRUD operations on plugins/components and bumps semantic versions.

**Complete documentation:** [README-auto-sync.md](./README-auto-sync.md)

**Key behaviors:**

- Returns exit code 0 (silent success, no commit interruption)
- Auto-stages modified manifest files
- Protects against double-bumping when commit fails and is retried
- Only operates on staged changes (`git diff --cached`)

### plugin_validator.py

Comprehensive validation tool for Claude Code plugins with token-based complexity measurement.

**Usage:**

```bash
# Validate single file
uv run plugins/plugin-creator/scripts/plugin_validator.py <path>

# Validate entire plugin
uv run plugins/plugin-creator/scripts/plugin_validator.py plugins/my-plugin

# Auto-fix issues
uv run plugins/plugin-creator/scripts/plugin_validator.py --fix <path>

# Validate only (no auto-fix)
uv run plugins/plugin-creator/scripts/plugin_validator.py --check <path>

# Verbose output
uv run plugins/plugin-creator/scripts/plugin_validator.py --verbose <path>

# CI mode (no color)
uv run plugins/plugin-creator/scripts/plugin_validator.py --no-color <path>
```

**Validates:**

- **Frontmatter:** YAML syntax, required fields, field types, tools/skills format
- **Plugin structure:** plugin.json schema, component paths, version consistency
- **Skill complexity:** Token-based metrics (4000 warning, 6400 error thresholds)
- **Internal links:** Markdown link validity, progressive disclosure
- **Completeness:** Required files, cross-references

**Auto-fixes:**

- YAML arrays → comma-separated strings
- Multiline descriptions → single-line quoted strings
- Unquoted descriptions with colons
- Removes `name:` field from plugin skills (Claude Code bug workaround)

**Error Codes:** 23 error codes across 9 validators - see ERROR_CODES.md

### create_plugin.py

Interactive plugin scaffolding tool.

**Usage:**

```bash
uv run plugins/plugin-creator/scripts/create_plugin.py
```

**Creates:**

- `.claude-plugin/` directory with `plugin.json`
- Optional `skills/`, `agents/`, `commands/` directories
- Self-validates before completion

### fix_tool_formats.py

Scans codebase for invalid tool format patterns in frontmatter and fixes them.

**Usage:**

```bash
# Scans ~/.claude and ~/repos recursively
uv run plugins/plugin-creator/scripts/fix_tool_formats.py
```

**Fixes:**

- YAML list → comma-separated string
- JSON array → comma-separated string

**Documentation:** [README.md](./README.md)

### count-skill-lines.sh

Counts lines in skills and identifies oversized ones (>500 lines warning, >800 lines critical).

**Usage:**

```bash
plugins/plugin-creator/scripts/count-skill-lines.sh <plugin-or-skill-path>
```

### validate-skill-structure.sh

Quality checks for skill structure beyond frontmatter validation.

**Validates:**

- Frontmatter presence and closure
- Name field format (lowercase, hyphens)
- Description length and trigger phrases
- Line count limits
- Progressive disclosure structure
- Internal link validity

**Exit codes:**

- 0: Pass
- 1: Fail (errors found)

### validate-task-file.sh

Validates refactoring task file format.

**Usage:**

```bash
plugins/plugin-creator/scripts/validate-task-file.sh <task-file-path>
```

## Pre-Commit Integration

Scripts integrated into `.pre-commit-config.yaml`:

| Hook ID                | Script                    | Trigger Pattern                                           | Purpose                                 |
| ---------------------- | ------------------------- | --------------------------------------------------------- | --------------------------------------- |
| `auto-sync-manifests`  | `auto_sync_manifests.py`  | `^plugins/`                                               | Auto-bump versions and update manifests |
| `plugin-validator`     | `plugin_validator.py`     | `^plugins/.*(SKILL\.md\|agents/.*\.md\|commands/.*\.md\|plugin\.json)$` | Comprehensive plugin validation with token metrics |

## Execution Requirements

All Python scripts require:

- Python 3.11+
- `uv` package manager
- No external dependencies (stdlib only)

Bash scripts require:

- Bash 5.1+
- Standard POSIX utilities

## Path Conventions

Scripts use absolute paths or `${CLAUDE_PLUGIN_ROOT}` variable when executed as commands.

When running manually:

- From repository root: `uv run plugins/plugin-creator/scripts/<script>`
- From scripts directory: `./<script>` (bash scripts only)

## Model Usage Guidelines

### When to Use These Scripts

**Use auto-sync-manifests behavior as reference when:**

- Documenting pre-commit hook behavior
- Explaining version bumping rules
- Troubleshooting manifest inconsistencies

**Use validation scripts when:**

- User reports frontmatter errors
- Creating or modifying plugin components
- Preparing plugins for marketplace submission

**Do NOT:**

- Execute scripts without explicit user request
- Modify scripts without reading complete source
- Assume behavior without verification

### Reading Documentation

Before referencing script behavior, the model MUST:

1. Read the script source code
2. Verify claims against actual implementation
3. Cite line numbers when describing behavior
4. Check associated README files for complete context

**Example citation:**

> "The script protects against double-bumping (lines 276-280 of auto_sync_manifests.py) by checking if plugin.json is already staged."

### Verification Protocol

When documenting script behavior:

1. Read the script file completely
2. Identify entry points and exit codes
3. Trace execution flow for the user's scenario
4. Cite specific line numbers for claims
5. Note any README documentation that exists

Do NOT:

- Assume behavior from script names
- Fabricate capabilities
- State opinions as facts
- Omit uncertainty when unable to verify
