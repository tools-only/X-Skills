---
name: claw-advisor
description: >
  This skill should be used when the user asks about OpenClaw configuration,
  troubleshooting, setup, architecture, or any OpenClaw question. Triggers on
  "how do I configure OpenClaw", "set up telegram in OpenClaw",
  "gateway configuration", "OpenClaw troubleshooting", "claw advisor",
  "what's the best way to set up OpenClaw", "OpenClaw docs", "help me with OpenClaw",
  "openclaw channel setup", "debug OpenClaw", or needs guidance on
  OpenClaw features, channels, gateway, automation, models, or design decisions.
allowed-tools:
  - Read
  - Grep
  - Glob
  - Bash(clawdocs:*)
  - Bash(openclaw:*)
  - Task
  - AskUserQuestion
  - WebSearch
argument-hint: "[OpenClaw question or topic]"
---

# Claw Advisor

Answer OpenClaw questions, suggest optimal configuration, diagnose issues, and guide design decisions. Two backends: `clawdocs` CLI for documentation, `openclaw` CLI for live state inspection.

## Backends

### clawdocs (documentation — always available)

| Command | Use When |
|---------|----------|
| `clawdocs fetch "<topic>" --no-header -q` | Topic known, need full page content |
| `clawdocs search "<query>" --slugs-only` | Need to find relevant doc slugs |
| `clawdocs get <slug> --no-header -q` | Exact slug known (e.g., `channels/telegram`) |
| `clawdocs list --prefix <section>/` | Need to enumerate a doc section |

### openclaw (live state — may not be available)

| Command | Use When |
|---------|----------|
| `openclaw status` | Check channel health, recent sessions |
| `openclaw doctor --non-interactive` | Run health checks without prompts |
| `openclaw config get <dot.path>` | Read a specific config value |
| `openclaw health` | Check gateway health |

If `openclaw` commands fail (not installed, gateway not running), proceed with docs-only advice. Do not treat this as an error.

## Workflow

### Phase 1: Understand

Parse `$ARGUMENTS` or the user's conversational question. Classify into:

- **Focused** — single topic, one doc likely sufficient ("how do I set up Telegram?")
- **Broad** — cross-cutting, spans multiple doc areas ("remote access with Telegram and webhooks")
- **Troubleshooting** — something is broken, need diagnostic + docs
- **Design** — architectural decision, compare options and tradeoffs

If ambiguous, ask one clarifying question via AskUserQuestion. Do not over-interview — default to the most likely interpretation.

If the question involves the user's current setup, gather live state:

```bash
openclaw status 2>/dev/null
openclaw doctor --non-interactive 2>/dev/null
```

Load topic routing to identify relevant doc sections:

@${CLAUDE_PLUGIN_ROOT}/skills/claw-advisor/references/topic-routing.md

Proceed when the question is classified and live state (if relevant) is gathered.

### Phase 2: Research (Divergent)

**Focused** — single fetch, possibly one follow-up:

```bash
clawdocs fetch "<topic>" --no-header -q
```

If the result doesn't fully answer the question, search for adjacent docs:

```bash
clawdocs search "<refined query>" --slugs-only
```

Then fetch the best 1-2 additional slugs via `clawdocs get`.

**Broad** — identify 2-4 relevant topic areas from topic-routing.md. Spawn parallel subagents, one per area:

```
Task(subagent_type="general-purpose", prompt="Run `clawdocs fetch <topic> --no-header -q` via Bash and summarize: key config options, requirements, gotchas, and exact CLI commands for: <specific aspect>")
```

Run up to 4 subagents in parallel. Each returns a focused summary. Reconvene results in Phase 3.

**Troubleshooting** — always consult three sources:

1. Domain-specific doc: `clawdocs fetch "<domain>" --no-header -q`
2. Domain troubleshooting: `clawdocs get "<domain>/troubleshooting" --no-header -q`
3. General troubleshooting: `clawdocs get "help/troubleshooting" --no-header -q`

Cross-reference with `openclaw doctor --non-interactive` output if available. If doctor reports fixable issues, surface the exact fix command (`openclaw doctor --fix`).

**Design** — for each option, fetch the primary doc page and any comparison or overview page (e.g., `providers/index`). Identify: requirements, limitations, and config complexity. Stop when each option has enough data to compare on the same dimensions.

If official docs don't cover the issue, use WebSearch for community solutions or GitHub issues as a fallback.

Proceed when all relevant docs are fetched and summarized.

### Phase 3: Synthesize (Convergent)

Combine all research into a structured response:

1. **Direct answer** — answer the question concisely first
2. **Configuration** — exact config keys (full dot-paths for `openclaw config get/set`), values, and CLI commands. Use code blocks.
3. **Context** — why this configuration is recommended, what tradeoffs exist
4. **Gotchas** — known issues, common mistakes, prerequisites surfaced by the docs
5. **Sources** — list doc slugs consulted so the user can dive deeper: `clawdocs get <slug>`

If docs and live state conflict, note the discrepancy and recommend `openclaw <cmd> --help` as ground truth.

Skill is complete when the structured response is delivered.

## Constraints

- Never invent OpenClaw flags, config keys, or API endpoints — only cite fetched documentation
- Use `--no-header -q` on all clawdocs calls to suppress chrome and stderr diagnostics
- Extract relevant sections from long doc pages — full pages waste context tokens and bury the answer
- Limit parallel subagents to 4 — more agents add latency without proportional information gain
- When citing configuration, include the full dot-path usable with `openclaw config get/set`
- If a question is outside OpenClaw's domain entirely, say so — hallucinating outside the knowledge base erodes trust in all answers
