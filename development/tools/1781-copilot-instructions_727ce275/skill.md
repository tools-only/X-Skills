# Copilot Coding Agent Instructions

Trust these instructions first. Only perform additional exploration if information is incomplete or found to be in error.

## Repository Overview

This repository is a **Claude Code Marketplace Plugin Collection** containing 22+ plugins that extend Claude's capabilities with specialized skills, commands, and agents. The plugins cover Python development, code quality, Git/CI-CD, AI/LLM tools, documentation, and agent orchestration.

**Project Type**: Plugin marketplace for Claude Code CLI
**Languages**: Markdown (primary for skills/commands/agents), Python 3.11+ (scripts)
**Package Manager**: uv (Astral's extremely fast Python package manager)
**Python Version**: 3.11+ required

## Directory Structure

```text
claude_skills/
├── .claude/                  # Claude Code local config (symlinked skills/commands/agents)
│   ├── agents/               # Agent markdown files
│   ├── commands/             # Command markdown files
│   └── skills/               # Skill directories (symlinks to plugins/)
├── .claude-plugin/           # Repository-level marketplace config
│   └── marketplace.json      # Plugin registry manifest
├── plugins/                  # All 22 plugins live here
│   └── {plugin-name}/
│       ├── .claude-plugin/
│       │   └── plugin.json   # Plugin manifest
│       ├── skills/           # Skill directories (SKILL.md + references/)
│       ├── commands/         # Command markdown files (optional)
│       ├── agents/           # Agent markdown files (optional)
│       └── README.md
├── sessions/                 # cc-sessions framework (task management)
├── scripts/                  # Utility scripts
├── marketplace.json          # Local marketplace config (in .claude-plugin/)
├── pyproject.toml            # Python project config (linting, type checking)
├── .pre-commit-config.yaml   # Pre-commit/prek hooks config
└── CLAUDE.md                 # AI-facing project instructions
```

## Build and Development Commands

### Environment Setup (Required First)

Always run these commands before any other operations:

```bash
# Install Python dependencies and create virtual environment
uv sync

# This creates .venv/ and installs all dependencies from pyproject.toml
# Must complete successfully before running any uv run commands
```

### Testing Plugins Locally

To test plugins during development:

```bash
# Option 1: Load specific plugins for this session
claude --plugin-dir ./plugins/python3-development --plugin-dir ./plugins/holistic-linting

# Option 2: Add local marketplace for persistent enable/disable
/plugin marketplace add ./.claude-plugin/marketplace.json

# Install plugins you need (--scope local keeps it gitignored)
/plugin install python3-development@jamie-bitflight-skills --scope local

# Disable when not needed
/plugin disable python3-development@jamie-bitflight-skills

# Re-enable when needed
/plugin enable python3-development@jamie-bitflight-skills
```

### Linting and Formatting

Always run linting before committing changes:

```bash
# Run ruff linter with auto-fix on specific files
uv run ruff check --fix path/to/file.py

# Run ruff formatter
uv run ruff format path/to/file.py

# Run ty type checking
uv run ty check path/to/file.py

# Run pre-commit hooks on specific files (preferred method)
uv run prek run --files path/to/file.py

# Run all pre-commit hooks on all files (slow, use sparingly)
uv run prek run --all-files
```

**Important Notes**:

- The repository uses `prek` (Rust-based pre-commit replacement), not `pre-commit`
- Some hooks may fail on first run due to network requirements (shellcheck-py)
- Run individual hooks when full pre-commit fails: `uv run prek run ruff --files <file>`

### Pre-commit Hook List

The following hooks run on commit:

| Hook                                       | Purpose                         |
| ------------------------------------------ | ------------------------------- |
| `ruff`                                     | Python linting with auto-fix    |
| `ruff-format`                              | Python formatting               |
| `markdownlint-cli2`                        | Markdown linting                |
| `prettier`                                 | YAML, JSON, Markdown formatting |
| `check-yaml`, `check-json`, `check-toml`   | Syntax validation               |
| `trailing-whitespace`, `end-of-file-fixer` | Whitespace normalization        |
| `ty-check`                                 | Python type checking            |
| `basedpyright`                             | Additional Python type checking |
| `conventional-pre-commit`                  | Commit message validation       |

## Plugin Structure

Each plugin must contain:

```text
plugins/{plugin-name}/
├── .claude-plugin/
│   └── plugin.json           # Required: plugin manifest
├── skills/
│   └── {skill-name}/
│       ├── SKILL.md          # Required: skill definition with YAML frontmatter
│       └── references/       # Optional: reference documentation
├── commands/                 # Optional: slash command definitions (.md files)
├── agents/                   # Optional: agent definitions (.md files)
└── README.md                 # Optional: user-facing documentation
```

### plugin.json Schema

```json
{
  "name": "plugin-name",
  "description": "Description with ACTION->TRIGGER->OUTCOME format",
  "version": "1.0.0",
  "skills": ["./skills/skill-name"],
  "commands": ["./commands"],
  "agents": ["./agents"]
}
```

### SKILL.md Format

Skills require YAML frontmatter:

```markdown
---
name: skill-name
description: Description with trigger conditions
---

# Skill Title

Content with rules, examples, and references...
```

## Coding Standards

### Python Files

- Target Python 3.11+
- Use `from __future__ import annotations` in all files
- Follow Google docstring convention
- Type hints required for all public functions
- Max line length: 120 characters
- Use native generics (list[str] not List[str])

### Markdown Files

- Skills are AI-facing documentation, NOT user documentation
- Use imperative language ("The model MUST...")
- Include XML tags for structured sections
- Cite all sources with URLs and access dates
- Use relative paths with `./` prefix for file references

### Commit Messages

Use Conventional Commits format:

```text
type(scope): description

Types: feat, fix, docs, style, refactor, test, chore
Scope: skill name or component
```

## Known Issues and Workarounds

### Pre-commit Network Failures

Some hooks require network access and may fail in restricted environments:

```bash
# If full pre-commit fails, run individual hooks
uv run prek run ruff ruff-format markdownlint-cli2 --files <changed-files>
```

### Ruff Format Trailing Comma

The project uses `skip-magic-trailing-comma = true`. Ruff may report files as needing reformatting for trailing comma differences. Run `uv run ruff format <file>` to fix.

### EXE003 Shebang Warning

Scripts using `uv run --script` shebang pattern trigger EXE003. This is intentionally ignored in `pyproject.toml`.

## Validation Checklist

Before submitting changes:

1. Run `uv sync` if dependencies changed
2. Run `uv run ruff check --fix` and `uv run ruff format` on Python files
3. Run `uv run ty check` on Python files
4. Verify plugin structure: `/plugin validate ./plugins/plugin-name`
5. Check markdown links are relative with `./` prefix
6. Ensure SKILL.md has valid YAML frontmatter
7. Use Conventional Commits format for commit messages

## File Locations Quick Reference

| Purpose              | Location                          |
| -------------------- | --------------------------------- |
| Linting config       | `pyproject.toml` [tool.ruff]      |
| Type checking config | `pyproject.toml` [tool.ty]        |
| Pre-commit hooks     | `.pre-commit-config.yaml`         |
| Markdown lint config | `.markdownlint-cli2.jsonc`        |
| Plugin registry      | `.claude-plugin/marketplace.json` |
| Project instructions | `CLAUDE.md`                       |
| Session framework    | `sessions/CLAUDE.sessions.md`     |

## Critical Reminders

- Skills are AI-facing documentation written for Claude, not humans
- Always use `uv run` prefix for Python commands
- Test new plugins with `claude --plugin-dir ./plugins/plugin-name`
- Paths in skills use `./` relative prefix
- CI runs via `.github/workflows/code-quality.yml` with quality gate pattern
- Environment setup for coding agent defined in `.github/workflows/copilot-setup-steps.yml`
