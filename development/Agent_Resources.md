---
name: Agent Resources
source: https://raw.githubusercontent.com/kasperjunge/agent-resources/main/AGENTS.md
original_path: AGENTS.md
source_repo: kasperjunge/agent-resources
category: development
subcategory: coding
tags: ['development']
collected_at: 2026-01-31T18:34:05.954675
file_hash: 263ca0123c9c49d2b8ffbdd032c61d093f5b4c90582c8a25b151bfc6cc497579
---

# Agent Resources

A package manager for AI agents. 

## Commands

This project uses `uv` for Python environment management. **Always use `uv run` to execute Python commands** to ensure they run in the correct virtual environment.

```bash
# Run tests
uv run pytest

# Run linters/formatters
uv run ruff check .
uv run ruff format .

# Run type checker
uv run ty check

# Test the CLI tools
uv run agr --help
uv run agrx --help
```

## Architecture
...

## agr.toml Format

The configuration file uses a flat array of dependencies:

```toml
dependencies = [
    {handle = "username/repo/skill", type = "skill"},
    {handle = "username/skill", type = "skill"},
    {path = "./local/skill", type = "skill"},
]
```

Each dependency has:
- `type`: Always "skill" for now
- `handle`: Remote GitHub reference (username/repo/skill or username/skill)
- `path`: Local path (alternative to handle)

Future: A `tools` section will configure which tools to sync to:
```toml
tools = ["claude", "cursor"]
```

## Code Style
...

## Boundaries

### Always Do
- agr and agrx should always be unified and synced.
- include in the plan to write tests for what is implemented
- Save all skills in `skills/` directory (not `.claude/skills/` which is gitignored)

### Ask First
...

### Never Do
...

## Security
...

# Docs

General
https://agentskills.io/
https://agents.md/

Claude Code:
https://code.claude.com/docs/en/skills
https://code.claude.com/docs/en/slash-commands
https://code.claude.com/docs/en/sub-agents
https://code.claude.com/docs/en/memory

Cursor:
https://cursor.com/docs/context/skills
https://cursor.com/docs/context/commands
https://cursor.com/docs/context/subagents
https://cursor.com/docs/context/rules

GitHub Copilot:
https://docs.github.com/en/copilot/concepts/agents/about-agent-skills

Codex:
https://developers.openai.com/codex/skills
https://developers.openai.com/codex/custom-prompts/

Open Code:
https://opencode.ai/docs/skills
https://opencode.ai/docs/commands/
https://opencode.ai/docs/agents/