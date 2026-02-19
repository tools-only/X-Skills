# TinyClaw - Multi-Agent Multi-Channel 24/7 AI Assistant Platform

**Research Date**: 2026-02-18
**GitHub**: <https://github.com/jlia0/tinyclaw>
**Version**: 0.0.5 (latest release 2026-02-15)
**License**: MIT
**Primary Language**: Shell (TypeScript runtime core)

---

## Overview

TinyClaw is an experimental multi-agent orchestration platform that runs named AI agents 24/7 across Discord, Telegram, and WhatsApp simultaneously. Agents are isolated processes (each running their own Claude Code or Codex CLI instance in a dedicated workspace), communicate peer-to-peer via bracket-tag mentions (`[@teammate: message]`), and are coordinated by a file-based queue system with atomic inbox/processing/outgoing directories. The system is designed for always-on personal AI assistants with team collaboration, live TUI visualization, and heartbeat-driven proactive task execution.

---

## Problem Addressed

| Problem | TinyClaw Solution |
|---------|-------------------|
| Claude Code runs in one terminal session -- no always-on availability | tmux-based daemon keeps agents running 24/7; `tinyclaw start/stop/restart` manages lifecycle |
| Single AI instance cannot specialize -- coding vs writing vs research require different prompts and models | Named isolated agents with per-agent provider, model, system prompt, and workspace directory |
| No cross-device access to AI assistant | Single logical assistant accessible from Discord (desktop), Telegram (mobile), WhatsApp (phone) simultaneously |
| Multi-agent orchestration requires central coordinator | Peer-to-peer team collaboration via `[@teammate: message]` tags -- no central orchestrator; agents directly pass work |
| AI sessions lose context across restarts | Persistent conversation history per agent maintained by the CLI; sessions survive `tinyclaw restart` |
| Race conditions in concurrent message handling | File-based queue with atomic `incoming/ → processing/ → outgoing/` state machine; no message corruption |
| Unknown senders can spam agents | Sender pairing allowlist with one-time code approval flow per channel |

---

## Key Statistics

| Metric | Value | Date Gathered |
|--------|-------|---------------|
| GitHub Stars | 2,124 | 2026-02-18 |
| GitHub Forks | 300 | 2026-02-18 |
| Open Issues | 54 | 2026-02-18 |
| Contributors | 4 | 2026-02-18 |
| Repository Size | ~17.6 MB | 2026-02-18 |
| Latest Version | 0.0.5 | 2026-02-18 |
| Created | 2026-02-09 | 2026-02-18 |
| Last Updated | 2026-02-18 | 2026-02-18 |
| License | MIT | 2026-02-18 |
| Primary Language | Shell + TypeScript | 2026-02-18 |
| Release Asset Size | 41.8 MB (bundle) | 2026-02-18 |
| Stability | Experimental | 2026-02-18 |
| Discord Community | Active (invite: discord.gg/jH6AcEChuD) | 2026-02-18 |
| Inspiration | OpenClaw by Peter Steinberger | 2026-02-18 |

**Note**: Created 2026-02-09, reached 2,124 stars in 9 days -- exceptionally fast community adoption for a pre-v1.0 experimental tool.

---

## Key Features

### Multi-Channel Architecture (3 Channels)

| Channel | Library | Notes |
|---------|---------|-------|
| Discord | `discord.js` v14 | Bot token setup, Message Content Intent required |
| Telegram | `node-telegram-bot-api` | BotFather token, standard bot |
| WhatsApp | `whatsapp-web.js` | QR code auth on first run; Linked Devices via phone |

All three channels share the same agent pool and conversation state. A user can start a conversation on Discord and continue on WhatsApp -- all messages route to the same persistent agent.

### Agent Isolation Model

Each agent is a fully isolated unit:

- **Dedicated workspace directory**: `~/tinyclaw-workspace/{agent_id}/`
- **Own AI CLI process**: Claude Code or Codex CLI launched per message
- **Independent conversation history**: Maintained by the CLI session
- **Custom configuration**: `.claude/` directory, `heartbeat.md`, `AGENTS.md`
- **Per-agent provider and model**: Mix Anthropic (Claude) and OpenAI (Codex) agents in the same deployment

Per-agent configuration in `.tinyclaw/settings.json`:

```json
{
  "agents": {
    "coder": {
      "name": "Code Assistant",
      "provider": "anthropic",
      "model": "sonnet",
      "working_directory": "/Users/me/tinyclaw-workspace/coder",
      "system_prompt": "You are a senior software engineer..."
    },
    "writer": {
      "name": "Technical Writer",
      "provider": "openai",
      "model": "gpt-5.3-codex",
      "working_directory": "/Users/me/tinyclaw-workspace/writer"
    }
  }
}
```

### Team Collaboration (Peer-to-Peer Agent Handoff)

Teams are named groups of agents with a designated leader. Work passes between agents via `[@teammate: message]` bracket tags in responses -- no central orchestrator involved.

**Team message flow:**

```text
User: "@dev fix the auth bug"
           │
           ▼
   Team: @dev
   Leader: @coder
           │
           ▼ "Fixed bug in auth.ts. [@reviewer: please review these changes]"
   @coder processes message
           │
           ▼ (new queue message generated)
   @reviewer: "Changes look good, approved!"
           │
           ▼ (all branches resolved)
   Aggregate → user receives both responses
```

**Supported message patterns** (documented in `docs/MESSAGE-PATTERNS.md`):

1. **Sequential handoff** -- single `[@teammate: message]` tag creates one hop
2. **Fan-out** -- multiple `[@teammate: message]` tags in one response, processed in parallel
3. **Backflow** -- agent mentions back the agent that originally mentioned it
4. **Cross-talk** -- agents after a fan-out mention each other before resolving
5. **Shared context** -- text outside bracket tags delivered to all mentioned agents

Team context is auto-detected even when messaging an agent directly -- if the agent belongs to a team, teammate mentions still trigger handoffs.

### File-Based Queue System

```text
~/.tinyclaw/queue/
├── incoming/     # new messages written atomically
├── processing/   # in-flight messages moved here
└── outgoing/     # completed responses staged here
```

Properties:
- **Atomic state transitions**: rename operations, not copies -- no partial reads
- **Parallel across agents**: different agents' messages processed concurrently
- **Sequential per agent**: preserves conversation order within each agent via per-agent promise chains
- **No external dependencies**: pure filesystem, no Redis, no database

### Live TUI Visualizer

React/Ink-based terminal UI (`tinyclaw team visualize [id]`):

- Agent cards with status (idle / active / done / error), provider, model, leader indicator
- Chain flow showing handoff arrows between agents
- Activity log of recent events with timestamps
- Status bar with queue depth and processing counts

Built with `ink` v6 (React-based TUI) and `ink-gradient`, `ink-spinner`.

### Heartbeat System

Agents can be configured with a `heartbeat.md` prompt file and an interval (seconds). The heartbeat fires proactively on schedule -- not in response to user messages. Default heartbeat prompt instructs the agent to check for pending tasks, errors, and unread messages, then take action if needed.

```json
{
  "monitoring": {
    "heartbeat_interval": 3600
  }
}
```

### Sender Pairing Allowlist

Unknown senders trigger a one-time pairing code flow:

1. First message from unknown sender generates a pairing code
2. Additional messages from same sender silently blocked while pending
3. User approves via `tinyclaw pairing approve <code>`
4. Subsequent messages from approved sender processed normally

### SOUL.md -- Agent Personality System

TinyClaw ships a `SOUL.md` template in each agent's workspace. This is a structured personality specification defining communication style, opinions, worldview, interests, and pet peeves. Agents with a filled-out `SOUL.md` produce responses consistent with that personality rather than defaulting to generic AI output.

---

## Technical Architecture

```text
┌─────────────────────────────────────────────────────────────┐
│                     Message Channels                         │
│         (Discord, Telegram, WhatsApp, Heartbeat)            │
│  discord-client.js  telegram-client.js  whatsapp-client.js  │
└────────────────────┬────────────────────────────────────────┘
                     │ Write message.json (atomic)
                     ↓
┌─────────────────────────────────────────────────────────────┐
│                   ~/.tinyclaw/queue/                         │
│                                                              │
│  incoming/          processing/         outgoing/           │
│  ├─ msg1.json  →   ├─ msg1.json   →   ├─ msg1.json        │
│  ├─ msg2.json       └─ msg2.json       └─ msg2.json        │
│  └─ msg3.json                                                │
└────────────────────┬────────────────────────────────────────┘
                     │ queue-processor.js
                     ↓
┌─────────────────────────────────────────────────────────────┐
│              Parallel Processing by Agent                    │
│                                                              │
│  Agent: coder              Agent: writer                    │
│  ┌─────────────────┐       ┌─────────────────┐             │
│  │ claude CLI      │       │ codex CLI        │             │
│  │ workspace/coder │       │ workspace/writer │             │
│  │ .claude/        │       │ AGENTS.md        │             │
│  │ heartbeat.md    │       │ heartbeat.md     │             │
│  │ AGENTS.md       │       │ SOUL.md          │             │
│  │ SOUL.md         │       └─────────────────┘             │
│  └─────────────────┘                                        │
│                                                              │
│  Peer-to-peer team handoff: [@teammate: message] tags       │
│  generate new messages in incoming/ without coordinator     │
└─────────────────────────────────────────────────────────────┘
                     │
                     ↓
┌─────────────────────────────────────────────────────────────┐
│              TUI Visualizer (team-visualizer.js)            │
│              React/Ink -- monitors ~/.tinyclaw/events/      │
└─────────────────────────────────────────────────────────────┘
```

### Technology Stack

| Component | Technology |
|-----------|-----------|
| Runtime | Node.js v18+ |
| Language | TypeScript 5.x (compiled to `dist/`) + Shell (orchestration scripts) |
| Discord client | discord.js v14 |
| Telegram client | node-telegram-bot-api |
| WhatsApp client | whatsapp-web.js v1.34 |
| TUI | React 19 + Ink v6 |
| JSON repair | jsonrepair v3 (auto-repairs settings.json with trailing commas/comments/BOM) |
| AI providers | Claude Code CLI (Anthropic) or Codex CLI (OpenAI) |
| Session management | tmux (daemon persistence) |

### Build System

```bash
npm run build         # tsc main + tsc visualizer
npm run queue         # node dist/queue-processor.js
npm run discord       # node dist/channels/discord-client.js
npm run visualize     # node dist/visualizer/team-visualizer.js
```

### Directory Structure

```text
tinyclaw/
├── .tinyclaw/            # Runtime data
│   ├── settings.json     # All configuration (agents, teams, channels)
│   ├── queue/            # File-based message queue
│   │   ├── incoming/
│   │   ├── processing/
│   │   └── outgoing/
│   ├── logs/             # Per-channel and per-component logs
│   ├── channels/         # Channel auth state (WhatsApp session)
│   ├── files/            # Uploaded files from channels
│   ├── pairing.json      # Sender allowlist state
│   ├── chats/            # Team conversation history (Markdown)
│   │   └── {team_id}/    # Timestamped .md files per conversation
│   └── events/           # Real-time event files for TUI visualizer
├── ~/tinyclaw-workspace/ # Agent workspaces (configurable path)
│   ├── coder/
│   │   ├── .claude/
│   │   ├── heartbeat.md
│   │   ├── AGENTS.md
│   │   └── SOUL.md
│   └── writer/
├── src/                  # TypeScript sources
├── dist/                 # Compiled output
├── lib/                  # Runtime shell scripts
└── scripts/              # Installation scripts
```

---

## Installation and Usage

### Prerequisites

- macOS, Linux, or Windows (WSL2)
- Node.js v18+
- tmux, jq
- Bash 4.0+ (macOS: `brew install bash`)
- Claude Code CLI (for Anthropic provider)
- Codex CLI (for OpenAI provider)

### Installation

```bash
# One-line install (recommended)
curl -fsSL https://raw.githubusercontent.com/jlia0/tinyclaw/main/scripts/remote-install.sh | bash

# From release bundle
wget https://github.com/jlia0/tinyclaw/releases/latest/download/tinyclaw-bundle.tar.gz
tar -xzf tinyclaw-bundle.tar.gz
cd tinyclaw && ./scripts/install.sh

# From source
git clone https://github.com/jlia0/tinyclaw.git
cd tinyclaw && npm install && ./scripts/install.sh
```

### Initial Setup

```bash
tinyclaw start    # Launches interactive setup wizard (7 steps)
```

Setup wizard covers: channel selection, bot tokens, workspace name, default agent config, AI provider, model, heartbeat interval.

### Core Management Commands

```bash
# Daemon lifecycle
tinyclaw start | stop | restart | status

# Agent management
tinyclaw agent list
tinyclaw agent add                             # interactive wizard
tinyclaw agent provider coder anthropic        # set provider
tinyclaw agent provider coder openai --model gpt-5.3-codex

# Team management
tinyclaw team list
tinyclaw team add                              # interactive wizard
tinyclaw team visualize dev                   # live TUI dashboard

# Message routing
tinyclaw send "Hello!"                        # default agent
tinyclaw send "@coder fix the auth bug"       # specific agent
```

### In-Chat Usage (Discord / Telegram / WhatsApp)

```text
@coder fix the authentication bug      → routes to coder agent
@dev fix the auth bug                  → routes to dev team leader
@writer document the API               → routes to writer agent
help me with this                      → routes to default agent

/agent                                 → list available agents
/team                                  → list available teams
@coder /reset                          → reset coder's conversation
```

### Team Configuration Example

```json
{
  "teams": {
    "dev": {
      "name": "Development Team",
      "agents": ["coder", "reviewer"],
      "leader_agent": "coder"
    }
  }
}
```

---

## Relevance to Claude Code Development

### Applications

1. **Peer-to-Peer Agent Communication Pattern**: TinyClaw's `[@teammate: message]` tag convention enables agents to spawn follow-on work without any central orchestrator. The queue processor merely detects these tags in responses and generates new queue entries. This decentralized handoff pattern eliminates the bottleneck and single point of failure of a coordinator agent, directly relevant to designing resilient multi-agent pipelines in Claude Code.

2. **File-Based Queue as Coordination Primitive**: The `incoming/ → processing/ → outgoing/` pattern achieves concurrent multi-agent coordination with zero infrastructure beyond a filesystem. Each agent's promise chain handles its own messages sequentially while different agents process in parallel -- a pattern immediately applicable to any Claude Code multi-agent orchestration needing reliable message delivery without Redis or a database.

3. **Isolated Workspace Per Agent**: Each TinyClaw agent has a separate directory with its own `.claude/` config, `CLAUDE.md`, and `AGENTS.md`. This workspace isolation model prevents context contamination between agents and allows per-agent customization of system prompts, memory files, and tools -- a pattern worth replicating when designing multi-agent Claude Code setups.

4. **SOUL.md as Personality Specification**: The structured SOUL.md template (vibe, worldview, opinions, interests, vocabulary) is a distinct approach to agent personality -- defining the AI's communication style declaratively, separate from task instructions. This separation of personality from capability is a pattern worth evaluating for agent design where consistent voice matters.

5. **tmux as Always-On Agent Substrate**: TinyClaw demonstrates that tmux is sufficient infrastructure for production 24/7 agent operation. No Kubernetes, no containers -- just named tmux sessions, cron for heartbeats, and shell scripts. This aligns with and extends the Claw Loop pattern already in this research collection.

### Patterns Worth Adopting

1. **Bracket Tag Handoff Convention**: `[@teammate: message]` is parseable with a simple regex, requires no changes to the AI provider, and works naturally within AI-generated prose. Any system that needs inter-agent communication via text output can adopt this convention.

2. **Atomic Queue State Machine**: The three-directory queue (incoming / processing / outgoing) uses filesystem rename operations for atomic state transitions -- no partial reads, no locking, no race conditions. This is a robust pattern for any file-based coordination system.

3. **Per-Channel, Per-Agent Log Separation**: `tinyclaw logs discord|telegram|whatsapp|queue|heartbeat|all` demonstrates the value of component-level log segregation for debugging multi-channel systems. Applying this to Claude Code multi-agent setups would significantly reduce debugging time.

4. **`jsonrepair` for Settings Tolerance**: Auto-repairing JSON configuration (trailing commas, comments, BOM characters) with a `.bak` backup before failing hard is a UX pattern that reduces user friction. Applicable wherever configuration files are human-edited.

5. **Team Chat History as Markdown Files**: Persisting team conversation chains as timestamped Markdown files in `~/.tinyclaw/chats/{team_id}/` creates an audit trail of agent-to-agent handoffs. This is more valuable than plain logs because the full context (user message, each agent's response, team metadata) is captured in a readable format.

### Integration Opportunities

1. **Claude Code Hook Integration**: TinyClaw's heartbeat system polls on a fixed interval. Claude Code hooks firing on specific events (task completion, context threshold, error conditions) could replace or augment this, enabling event-driven agent activation rather than time-based polling.

2. **SOUL.md as Skill Extension Point**: The SOUL.md personality specification could be formalized as a Claude Code skill that auto-loads into an agent's CLAUDE.md during setup, making agent personality configuration a first-class skill alongside task-focused skills.

3. **Queue Visualizer as MCP Tool**: The TUI visualizer reads `~/.tinyclaw/events/` -- this event bus architecture could be exposed as an MCP server, allowing external tools to observe and interact with running agent teams.

### Considerations

1. **Experimental Stability**: The project is explicitly marked experimental (v0.0.5, created 2026-02-09). The 54 open issues and rapid release cadence indicate active stabilization. Production use carries risk of breaking changes between releases.

2. **Shell + TypeScript Split**: The primary GitHub language is Shell, but the runtime core is TypeScript. This dual-language architecture can complicate contributions and debugging -- shell scripts orchestrate TypeScript processes, and bugs can exist in either layer.

3. **Claude Code CLI Dependency**: TinyClaw does not use the Claude API directly -- it invokes the Claude Code CLI (or Codex CLI) as a subprocess. This means it inherits Claude Code's session management, rate limits, and authentication model, and is not compatible with raw API access.

4. **WhatsApp ToS Risk**: `whatsapp-web.js` automates the WhatsApp Web interface, which is an unofficial API. WhatsApp has historically restricted unofficial automation clients. The README explicitly notes using "existing subscriptions without breaking ToS" for AI providers but makes no similar claim for WhatsApp itself.

5. **No Test Suite**: No test directory or test scripts are present in the repository. For a multi-channel daemon managing persistent AI agent sessions, the lack of automated tests makes regression risk high during the rapid development phase.

---

## References

1. **TinyClaw GitHub Repository** - <https://github.com/jlia0/tinyclaw> (accessed 2026-02-18)
2. **TinyClaw README.md** - <https://raw.githubusercontent.com/jlia0/tinyclaw/main/README.md> (accessed 2026-02-18)
3. **TinyClaw TEAMS.md** - <https://raw.githubusercontent.com/jlia0/tinyclaw/main/docs/TEAMS.md> (accessed 2026-02-18)
4. **TinyClaw AGENTS.md** - <https://raw.githubusercontent.com/jlia0/tinyclaw/main/docs/AGENTS.md> (accessed 2026-02-18)
5. **TinyClaw package.json** - <https://raw.githubusercontent.com/jlia0/tinyclaw/main/package.json> (accessed 2026-02-18)
6. **GitHub API Repository Metadata** - `https://api.github.com/repos/jlia0/tinyclaw` -- stars: 2124, forks: 300, language: Shell, license: MIT (accessed 2026-02-18)
7. **GitHub API Latest Release** - v0.0.5, published 2026-02-15, asset: tinyclaw-bundle.tar.gz (41.8 MB), release note: "Add `agent provider` command" (accessed 2026-02-18)
8. **GitHub API Contributors** - 4 contributors (Link header page count, accessed 2026-02-18)
9. **OpenClaw by Peter Steinberger** - Credited as inspiration in TinyClaw README (accessed 2026-02-18 via [1])

---

## Freshness Tracking

| Field | Value |
|-------|-------|
| Version Documented | 0.0.5 |
| Repository Created | 2026-02-09 |
| Latest Release | 2026-02-15 |
| Stars at Research Date | 2,124 |
| Open Issues | 54 |
| Stability | Experimental |
| Research Date | 2026-02-18 |
| Next Review | 2026-05-18 |

### Update Triggers

- Version bump to v0.1.0 or v1.0.0 (stability milestone)
- Addition of new channel integrations beyond Discord/Telegram/WhatsApp
- Addition of new AI provider support beyond Anthropic and OpenAI
- Introduction of test suite (significant quality signal for experimental project)
- Breaking changes to SOUL.md format, team configuration schema, or queue protocol
- OpenClaw integration or collaboration that affects architecture
- Community Discord activity indicating direction change or feature freeze
- Stars crossing 5K or 10K (adoption signal for risk assessment update)
