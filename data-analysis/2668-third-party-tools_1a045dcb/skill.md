---
title: "Third-Party Tools for Claude Code"
description: "Community tools for token tracking, session management, configuration, and alternative UIs"
tags: [reference, integration, plugin]
---

# Third-Party Tools for Claude Code

> Community tools for token tracking, session management, configuration, and alternative UIs.
>
> **Last verified**: February 2026

## Table of Contents

1. [About This Page](#about-this-page)
2. [Token & Cost Tracking](#token--cost-tracking)
3. [Session Management](#session-management)
4. [Configuration Management](#configuration-management)
5. [Alternative UIs](#alternative-uis)
6. [Multi-Agent Orchestration](#multi-agent-orchestration)
7. [Plugin Ecosystem](#plugin-ecosystem)
8. [Known Gaps](#known-gaps)
9. [Recommendations by Persona](#recommendations-by-persona)

---

## About This Page

This page catalogs **community-built tools that extend Claude Code**. Each tool has been verified against its public repository or package registry. Only tools with a public source (GitHub, npm, PyPI) are included.

**What this page is NOT**:
- Not a list of AI tools that complement Claude Code (see [AI Ecosystem](./ai-ecosystem.md))
- Not DIY monitoring scripts (see [Observability](./observability.md))
- Not MCP server recommendations (see [MCP Servers Ecosystem](./mcp-servers-ecosystem.md))

---

## Token & Cost Tracking

### ccusage

The most mature cost tracking tool for Claude Code. Parses local session data to produce cost reports by day, month, session, or 5-hour billing window.

| Attribute | Details |
|-----------|---------|
| **Source** | [npm: ccusage](https://www.npmjs.com/package/ccusage) / [ccusage.com](https://ccusage.com) |
| **Install** | `bunx ccusage` (fastest) or `npx ccusage` |
| **Language** | TypeScript (Node.js 18+) |
| **Version** | 18.x (actively maintained) |

**Key features**:

- `ccusage daily` / `ccusage monthly` / `ccusage session` - aggregated cost reports
- `ccusage blocks --live` - real-time monitoring against 5-hour billing windows
- `--breakdown` flag for per-model cost split (Opus/Sonnet/Haiku)
- `--since` / `--until` date filtering
- JSON output (`--json`) for programmatic access
- Offline mode with cached pricing data
- MCP server integration (`@ccusage/mcp`)
- macOS widget (`ccusage-widget`) and [Raycast extension](https://www.raycast.com/nyatinte/ccusage)

**Limitations**: Relies on local JSONL parsing; cost estimates may differ from official Anthropic billing. No team aggregation without manual log merging.

> **Cross-ref**: The main guide covers basic ccusage commands at [ultimate-guide.md Section 2.4](./ultimate-guide.md) (cost monitoring).
> For DIY cost tracking with hooks, see [Observability](./observability.md).

---

### ccburn

A Python TUI for visual token burn-rate tracking. Displays charts showing consumption rate relative to Claude's billing windows.

| Attribute | Details |
|-----------|---------|
| **Source** | [GitHub: JuanjoFuchs/ccburn](https://github.com/JuanjoFuchs/ccburn) / [Blog post](https://juanjofuchs.github.io/ai-development/2026/01/13/introducing-ccburn-visual-token-tracking.html) |
| **Install** | `pip install ccburn` |
| **Language** | Python 3.10+ (Rich + Plotext) |

**Key features**:

- Terminal charts showing token consumption over time
- Burn-rate indicators (on-track / slow-down warnings)
- Compact display mode
- Visual budget tracking against limits

**Limitations**: Python-only ecosystem. Smaller community than ccusage. No MCP integration.

**When to choose ccburn over ccusage**: If you prefer visual burn-rate charts over tabular reports, or if your toolchain is Python-based.

---

### RTK (Rust Token Killer)

A CLI proxy that filters command outputs **before** they reach Claude's context. 446 stars, 38 forks, 700+ upvotes on r/ClaudeAI.

| Attribute | Details |
|-----------|---------|
| **Source** | [GitHub: rtk-ai/rtk](https://github.com/rtk-ai/rtk) |
| **Website** | [rtk-ai.app](https://www.rtk-ai.app/) |
| **Install** | `brew install rtk-ai/tap/rtk` or `cargo install rtk` |
| **Language** | Rust (standalone binary) |
| **Version** | v0.16.0 |

**Key features**:

- `rtk git log` (92% reduction), `rtk git status` (76% reduction), `rtk git diff` (56% reduction)
- `rtk vitest run`, `rtk prisma`, `rtk pnpm` (70-90% reduction)
- `rtk python pytest`, `rtk go test` (multi-language support)
- `rtk cargo test/build/clippy` (Rust toolchain)
- `rtk init` - hook-first install, `rtk tree` - project structure, `rtk learn` - interactive learning
- `rtk gain` - token savings analytics (SQLite tracking)
- `rtk discover` - find missed optimization opportunities

**When to choose RTK vs ccusage/ccburn**:

- RTK **reduces** token consumption (preprocessing)
- ccusage/ccburn **monitor** it (postprocessing)
- Use both together for maximum efficiency

**Limitations**: Rapid development cadence (30 releases in 23 days). Not suitable for interactive commands or very small outputs.

> **Cross-ref**: Full docs at [ultimate-guide.md Section 9](./ultimate-guide.md#command-output-optimization-with-rtk)

---

## Session Management

### claude-code-viewer

A web-based UI for browsing and reading Claude Code conversation history (JSONL files).

| Attribute | Details |
|-----------|---------|
| **Source** | [GitHub: d-kimuson/claude-code-viewer](https://github.com/d-kimuson/claude-code-viewer) / [npm: @kimuson/claude-code-viewer](https://www.npmjs.com/package/@kimuson/claude-code-viewer) |
| **Install** | `npx @kimuson/claude-code-viewer` or `npm install -g @kimuson/claude-code-viewer` |
| **Language** | TypeScript (Node.js 18+) |
| **Version** | 0.5.x |

**Key features**:

- Project browser with session counts and metadata
- Full conversation display with syntax highlighting
- Tool usage results inline
- Real-time updates via Server-Sent Events (auto-refreshes when files change)
- Responsive design (desktop + mobile)

**Limitations**: Read-only (cannot edit or resume sessions). No cost data. Requires existing `~/.claude/projects/` history.

> **Cross-ref**: For session search from the CLI, see [session-search.sh](../examples/scripts/session-search.sh) in [Observability](./observability.md).

---

ti### Entire CLI

Agent-native platform for Git-integrated session capture with rewindable checkpoints and governance layer.

| Attribute | Details |
|-----------|---------|
| **Source** | [GitHub: entireio/cli](https://github.com/entireio/cli) / [entire.io](https://entire.io) |
| **Install** | See GitHub (platform launched Feb 2026, early access) |
| **Language** | TypeScript |
| **Founded** | February 2026 by Thomas Dohmke (ex-GitHub CEO), $60M funding |

**Key features:**

- **Session Capture**: Automatic recording of AI agent sessions (Claude Code, Gemini CLI) with full context
- **Rewindable Checkpoints**: Restore to any session state with prompts + reasoning + file changes
- **Governance Layer**: Permission system, human approval gates, audit trails for compliance
- **Agent Handoffs**: Preserve context when switching between agents (Claude → Gemini)
- **Git Integration**: Stores checkpoints on separate `entire/checkpoints/v1` branch (no history pollution)
- **Multi-Agent Support**: Works with multiple AI agents simultaneously with context sharing

**Use cases:**

| Scenario | Why Entire CLI |
|----------|---------------|
| **Compliance (SOC2, HIPAA)** | Full audit trail: prompts → reasoning → outputs |
| **Multi-agent workflows** | Context preserved across agent switches |
| **Debugging AI decisions** | Rewind to checkpoint, inspect reasoning |
| **Governance** | Approval gates before production changes |
| **Team handoffs** | Resume sessions with full context |

**vs claude-code-viewer:**

| Feature | claude-code-viewer | Entire CLI |
|---------|-------------------|-----------|
| **Purpose** | Read-only history viewing | Active session management + replay |
| **Replay** | No | Yes (rewind to checkpoints) |
| **Context** | Conversation only | Prompts + reasoning + file states |
| **Governance** | No | Yes (approval gates, permissions) |
| **Multi-agent** | No | Yes (agent handoffs) |
| **Overhead** | None | ~5-10% storage |

**When to choose Entire over claude-code-viewer:**

- ✅ Need session replay/rewind functionality
- ✅ Enterprise compliance requirements (audit trails)
- ✅ Multi-agent workflows (Claude + Gemini)
- ✅ Governance gates (approval before deploy)
- ❌ Just want to browse history → Use claude-code-viewer (lighter)

**Limitations:**

- Very new (launched Feb 10-12, 2026) - limited production feedback
- Enterprise-focused (may be complex for solo developers)
- Storage overhead (~5-10% of project size for session data)
- macOS/Linux only (Windows via WSL)
- Early stage (v1.x) - expect API changes

> **Cross-ref**: Full Entire workflow with examples at [AI Traceability Guide](./ai-traceability.md#51-entire-cli). For compliance use cases, see [Security Hardening](./security-hardening.md).

---

## Configuration Management

### claude-code-config

A TUI for managing `~/.claude.json` configuration, focused on MCP server management.

| Attribute | Details |
|-----------|---------|
| **Source** | [GitHub: joeyism/claude-code-config](https://github.com/joeyism/claude-code-config) |
| **Install** | `pip install claude-code-config` |
| **Language** | Python (Textual TUI) |

**Key features**:

- Visual MCP server management (add, edit, remove)
- Configuration file editing with validation
- TUI navigation for `~/.claude.json` structure

**Limitations**: Limited to `~/.claude.json` scope. Does not manage `.claude/settings.json`, hooks, or slash commands.

---

### AIBlueprint

A CLI that scaffolds pre-configured Claude Code setups with hooks, commands, statusline, and workflow automation.

| Attribute | Details |
|-----------|---------|
| **Source** | [GitHub: Melvynx/aiblueprint](https://github.com/Melvynx/aiblueprint) |
| **Install** | `npx aiblueprint-cli` |
| **Language** | TypeScript |

**Key features**:

- Pre-built security hooks
- Custom command templates
- Statusline configuration
- Workflow automation presets

**Limitations**: Opinionated configuration choices. Some features require a premium tier. Does not read existing config (scaffolds from scratch).

> **Cross-ref**: For manual Claude Code configuration, see [ultimate-guide.md Section 4](./ultimate-guide.md) (CLAUDE.md, settings, hooks, commands).

---

## Alternative UIs

### Claude Chic

A styled terminal UI for Claude Code built on Anthropic's claude-agent-sdk. Replaces the default Claude Code TUI with a visually enhanced experience.

| Attribute | Details |
|-----------|---------|
| **Source** | [Blog: matthewrocklin.com](https://matthewrocklin.com/introducing-claude-chic/) / [PyPI: claudechic](https://pypi.org/project/claudechic/) |
| **Install** | `uvx claudechic` |
| **Language** | Python (Textual + claude-agent-sdk) |
| **Status** | Alpha |

**Key features**:

- Color-coded messages (orange: user, blue: Claude, grey: tools)
- Collapsible tool usage blocks
- Git worktree management from within the UI
- Multiple agents in a single window
- `/diff` viewer, vim keybindings (`/vim`), shell commands (`!ls`)
- Proper Markdown rendering with streaming

**Limitations**: Alpha status - expect breaking changes. Python dependency chain. Requires claude-agent-sdk. macOS/Linux only.

---

### Toad

A universal terminal frontend for AI coding agents. Supports Claude Code alongside Gemini CLI, OpenHands, Codex, and 12+ other agents via the Agent Client Protocol (ACP).

| Attribute | Details |
|-----------|---------|
| **Source** | [GitHub: batrachianai/toad](https://github.com/batrachianai/toad) / [willmcgugan.github.io/toad-released](https://willmcgugan.github.io/toad-released/) |
| **Install** | `curl -fsSL batrachian.ai/install \| sh` or `uv tool install -U batrachian-toad --python 3.14` |
| **Author** | Will McGugan (creator of Rich & Textual) |
| **Language** | Python (Textual) |

**Key features**:

- Unified interface across 12+ agent CLIs
- Full shell integration with tab completion
- `@` file context injection with fuzzy search
- Side-by-side diffs with syntax highlighting
- Jupyter-inspired block navigation
- Flicker-free character-level rendering

**Limitations**: macOS/Linux only (Windows via WSL). Agent support varies by ACP compatibility. No built-in session persistence yet (on roadmap).

---

### Conductor

A macOS desktop app for orchestrating multiple Claude Code instances in parallel using git worktrees.

| Attribute | Details |
|-----------|---------|
| **Source** | [conductor.build](https://docs.conductor.build) |
| **Install** | Download from [conductor.build](https://docs.conductor.build) |
| **Platform** | macOS only (Windows/Linux planned) |

**Key features**:

- Parallel Claude Code agents on separate git worktrees
- Unified diff viewer for all agent changes
- GitHub and Linear integration (inject issues as context)
- MCP support and slash commands (`/research`)
- Planning mode with bulk context injection

**Limitations**: macOS only (as of Feb 2026). Proprietary (not open source). Overlaps with multi-agent orchestration (see below).

---

### Claude Code GUI (VS Code Extension)

A third-party VS Code extension (not Anthropic's official extension) that adds a graphical layer on top of Claude Code.

| Attribute | Details |
|-----------|---------|
| **Source** | [VS Code Marketplace: MaheshKok.claude-code-gui](https://marketplace.visualstudio.com/items?itemName=MaheshKok.claude-code-gui) |
| **Install** | VS Code Marketplace → search "Claude Code GUI" |

**Note**: This is **not** the official [Claude Code for VS Code](https://marketplace.visualstudio.com/items?itemName=anthropic.claude-code) extension by Anthropic. The official extension provides inline diffs, @-mentions, and plan review directly in the editor.

**Limitations**: Third-party, not Anthropic-maintained. Feature set may overlap with or lag behind the official extension.

---

## Multi-Agent Orchestration

This section covers tools for running **multiple Claude Code instances in parallel**. For detailed documentation, see:

- **[AI Ecosystem](./ai-ecosystem.md)** - Gas Town, multiclaude, agent-chat, claude-squad
- **[Ultimate Guide Section 9](./ultimate-guide.md)** - Multi-instance workflows, git worktrees, orchestration frameworks

**Quick reference**:

| Tool | Type | Key Feature |
|------|------|-------------|
| [Gas Town](https://github.com/steveyegge/gastown) | Multi-agent workspace | Steve Yegge's agent-first workspace manager |
| [multiclaude](https://github.com/dlorenc/multiclaude) | Multi-agent spawner | tmux + git worktrees (383+ stars) |
| [agent-chat](https://github.com/justinabrahms/agent-chat) | Monitoring UI | Real-time SSE monitoring for Gas Town/multiclaude |
| [Conductor](#conductor) | Desktop app | macOS parallel agents (also listed above) |

---

## Plugin Ecosystem

Claude Code's plugin system supports community-built extensions. For detailed documentation:

- **[Ultimate Guide Section 8](./ultimate-guide.md)** - Plugin system, commands, installation
- **[claude-plugins.dev](https://claude-plugins.dev)** - 11,989 plugins, 63,065 skills indexed
- **[claudemarketplaces.com](https://claudemarketplaces.com)** - Auto-scan GitHub for marketplace plugins
- **[agentskills.io](https://agentskills.io)** - Open standard for agent skills (26+ platforms)

---

## Known Gaps

As of February 2026, the community tooling ecosystem has notable gaps:

| Gap | Description |
|-----|-------------|
| **Visual skills editor** | No GUI for creating/editing `.claude/skills/` — must edit YAML/Markdown manually |
| **Visual hooks editor** | No GUI for managing hooks in `settings.json` — requires JSON editing |
| **Unified admin panel** | No single dashboard combining config, sessions, cost, and MCP management |
| **Session replay** | ✅ **FILLED**: Entire CLI (launched Feb 2026) provides rewindable checkpoints with full context replay |
| **Agent-native issue tracking** | No established tool for markdown-based, git-committable issue tracking with Claude Code. [fp.dev](https://fp.dev/) is an early-stage solution (local-first, `/fp-plan` + `/fp-implement` skills, diff viewer) but lacks adoption signals and requires Apple Silicon for the desktop app. The Tasks API covers state persistence but issues aren't git-committable. |
| **Per-MCP-server profiler** | No way to measure token cost attributable to each MCP server individually |
| **Cross-platform config sync** | No tool syncs Claude Code config across machines (must manual copy `~/.claude/`) |

---

## Recommendations by Persona

| Persona | Recommended Tools | Rationale |
|---------|-------------------|-----------|
| **Solo developer** | ccusage + claude-code-viewer | Cost awareness + session history review |
| **Small team (2-5)** | ccusage + Conductor or multiclaude | Cost tracking + parallel development |
| **Enterprise** | ccusage (MCP) + custom dashboards | Programmatic cost data + audit trails |
| **Python-centric** | ccburn + Claude Chic | Native Python ecosystem tools |
| **Multi-agent user** | Toad or Conductor | Unified agent management |
| **Config-heavy setup** | claude-code-config + AIBlueprint | TUI config management + scaffolding |

---

## Related Resources

- [Observability](./observability.md) - DIY session monitoring, logging hooks, cost tracking scripts
- [AI Ecosystem](./ai-ecosystem.md) - Complementary AI tools (Perplexity, Gemini, NotebookLM)
- [MCP Servers Ecosystem](./mcp-servers-ecosystem.md) - Validated community MCP servers
- [Architecture](./architecture.md) - How Claude Code works internally
- [Ultimate Guide Section 8](./ultimate-guide.md) - Plugin system and marketplaces
