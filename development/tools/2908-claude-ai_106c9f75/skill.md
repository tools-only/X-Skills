# Claude AI (claude.ai / Desktop App)

> Last verified: 2026-02-17

## Models

Default model varies by plan: Opus 4.6 for Pro/Max/Team, Sonnet 4.5 for Free. See [models.md](models.md) for full list.

## Plans

| Tier | Key Limits |
|------|-----------|
| Free | File creation, connectors, skills, longer conversations (expanded Jan 2026) |
| Pro | Higher usage limits, all features |
| Max | Highest limits, Cowork access |
| Team | Collaboration features, admin controls |
| Enterprise | SSO, audit logs, custom retention |

## Skills

Upload via Settings > Capabilities > Skills. Format: zip archive containing SKILL.md + optional scripts/, references/, assets/. Skills trigger based on YAML frontmatter `description` field. Available on all plans including Free.

## MCP Connectors

Connected via UI (Settings > Connectors). Current first-party connectors include Gmail, Google Chat, Google Drive, Google Docs, Google Sheets, Google Tasks, Cloudflare, Mermaid Chart. Third-party connectors available through partner integrations.

## Interactive Apps

Launched Jan 26, 2026. Powered by "MCP Apps" protocol extension. Nine launch partners:

- Slack, Canva, Figma, Box, Clay, Asana, Amplitude, Hex, Monday.com

These deliver interactive UI within the Claude conversation. Different from MCP connectors — they're rich, embedded app experiences. Available on Pro, Max, Team, Enterprise.

## Artifacts

Render inline in conversation:

| Type | Description |
|------|-------------|
| HTML | Full web pages with CSS/JS |
| React/JSX | Interactive components |
| Mermaid | Diagrams and charts |
| SVG | Vector graphics |
| Markdown | Formatted documents |
| Code | Syntax-highlighted, copyable |

## File Creation

Can create and download: DOCX, PPTX, XLSX, PDF, code files, images (via tool use). Output path: `/mnt/user-data/outputs/`. Use `present_files` tool to share files with the user.

## Memory

Automatic cross-conversation memory from chat history. Claude remembers facts, preferences, and context from previous conversations. No explicit memory management by user — it just works.

## Search

- **Web search**: Built-in with citations. Automatic when current information is needed.
- **Conversation search**: `conversation_search` and `recent_chats` internal tools for searching past conversations. Not accessible from external tools or APIs.

## Beta Features

- **Claude in Chrome**: Browser automation agent. Controls Chrome tabs, fills forms, navigates. Beta access.
- **Claude in Excel**: Spreadsheet automation agent. Beta access.

## Container Environment

Ubuntu-based sandbox. Python available for skill script execution. Key paths:

- `/mnt/user-data/outputs/` — file output directory
- `/mnt/skills/` — loaded skill files

No persistent bash shell. No local filesystem access beyond the container. Cannot install arbitrary packages. No git, no deployment tools, no CLI access to external services.

## What Claude AI Cannot Do

- Access the user's local filesystem
- Run bash commands on the user's machine
- Deploy code (no wrangler, vercel, etc.)
- Use subagents or Agent Teams
- Access Claude Code conversation logs
- Run persistent background processes
