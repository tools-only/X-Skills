# Claude Code (Terminal / IDE)

> Last verified: 2026-02-17

## Overview

CLI tool for software engineering. Runs in the terminal with full local filesystem and bash access. Current version: ~v2.1.42 (Feb 13, 2026).

## Plans

Available on all Claude plans (Pro, Max, Team, Enterprise). Some features require Max plan (see Agent Teams, Cowork below).

## Skills

Directory-based, no packaging needed for local use:

- Global: `~/.claude/skills/<skill-name>/SKILL.md`
- Project: `.claude/skills/<skill-name>/SKILL.md`
- Install: `/plugin install <name>` or `/plugin marketplace add <repo>`

Skills auto-create `/skill-name` slash commands.

## MCP

Configure via `claude mcp add <name>` or edit settings.json directly. Supports stdio and SSE transports. Tools appear as `mcp__<server>__<tool>` format. MCP tools are inherited by subagents when the `tools` field is omitted from agent config.

## Agent Teams (Experimental)

Launched Feb 5, 2026. Multi-instance orchestration where multiple Claude Code sessions coordinate as a team.

- **Shared task lists**: TaskCreate, TaskUpdate, TaskList, TaskGet tools
- **Inter-agent messaging**: Teammates communicate directly
- **Display**: In-process (Shift+Up/Down) or split panes (tmux/iTerm2)
- **Enable**: `CLAUDE_CODE_EXPERIMENTAL_AGENT_TEAMS=1`
- **Requires**: Max plan
- **Use cases**: Research + review in parallel, cross-layer coordination, parallel debugging hypotheses

Each teammate gets its own context window. Token usage scales with teammate count.

## Cowork

Launched Jan 2026. Claude Code for non-developers. Provides local folder access with file watching.

- macOS only
- Requires Max plan
- Designed for knowledge workers, not just software engineers

## Subagents

Fork work to specialised agents via the Task tool:

| Built-in | Purpose |
|----------|---------|
| Explore | Read-only codebase exploration |
| Plan | Architecture planning without edits |
| general-purpose | Flexible multi-tool agent |

Custom agents defined in `.claude/agents/` directory. Subagents inherit MCP tools when the `tools` field is omitted from the agent config.

## Hooks

Event-driven automation configured in settings or `.claude/` directory:

| Event | Fires when |
|-------|-----------|
| PreToolUse | Before a tool is called |
| PostToolUse | After a tool completes |
| Notification | On notifications |
| TeammateIdle | Agent Teams: teammate finishes work |
| TaskCompleted | Agent Teams: task marked complete |

## Memory System

Hierarchical CLAUDE.md files loaded from current directory + all parent directories:

- `~/CLAUDE.md` — global context (always loaded)
- `~/project/CLAUDE.md` — project context
- `~/project/src/CLAUDE.md` — directory-specific context

Additional memory:
- `~/.claude/rules/` — global rules (always loaded)
- `.claude/rules/` — project rules
- Memory frontmatter scopes: `user`, `project`, `local`

## Available Tools

| Tool | Purpose |
|------|---------|
| Read, Write, Edit | File operations |
| Glob, Grep | Search files and content |
| Bash | Shell command execution |
| WebFetch, WebSearch | Web access |
| Task | Delegate to subagents |
| ToolSearch | Discover deferred/MCP tools |
| AskUserQuestion | Get user input |

Plus any MCP tools from connected servers.

## What Claude Code Cannot Do

- Render artifacts (HTML, React) inline
- Access Claude AI conversation history (no API exists)
- Use Interactive Apps (Slack, Canva, etc.)
- Run without a terminal/IDE environment
- Create DOCX/PPTX/XLSX natively (needs libraries)
- Search past Claude AI conversations
