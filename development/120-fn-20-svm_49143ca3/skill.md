# fn-20-svm Codex Multi-agents Parity

## Status: BLOCKED

Blocked on native custom subagent support in Codex.

## Overview

Codex 0.88.0 (Jan 21, 2026) added experimental Multi-agents feature. This epic tracks bringing flow-next's subagent capabilities to Codex once the platform supports custom agent definitions.

## Current Codex Multi-agents

Enable via `/experiments` → `[x] Multi-agents`

### Collab tools (orchestrator only)
```
spawn_agent(prompt, agent_type)  → spawn worker
send_input(thread_id, message)   → send follow-ups
send_input(..., interrupt=true)  → stop + redirect
wait(thread_ids)                 → block until done
close_agent(thread_id)           → terminate
```

### Built-in agent_type options
| Type | Purpose |
|------|---------|
| `default` | Inherits parent config |
| `orchestrator` | Coordination-only, spawns workers |
| `worker` | Task-executing agent |

## Gap Analysis

| Claude Code | Codex | Gap |
|-------------|-------|-----|
| 12+ custom agents | 3 generic roles | No custom definitions |
| Per-agent tools/model | Same config for all | No per-agent config |
| `agents/*.md` definitions | Just role type | No agent storage |
| Task tool + subagent_type | spawn_agent + agent_type | Limited types |

### flow-next agents not portable to Codex
- `practice-scout` - best practices research
- `docs-scout` - documentation lookup
- `repo-scout` - codebase patterns
- `github-scout` - GitHub code search
- `context-scout` - RepoPrompt exploration
- `docs-gap-scout` - docs update detection
- `memory-scout` - memory system search
- `epic-scout` - epic dependencies
- `flow-gap-analyst` - requirements analysis
- `quality-auditor` - code review
- `plan-sync` - spec synchronization
- `worker` - task implementation

## Workaround (not recommended)

The [codex-subagents-mcp](https://github.com/leonardsellem/codex-subagents-mcp) pattern:
1. Create temp workdir with custom AGENTS.md (scout persona)
2. Spawn worker in that dir via `codex exec --profile`
3. Worker inherits persona from AGENTS.md

Problems:
- Fragile temp directory management
- No native tool restrictions per agent
- No model override per agent
- Extra overhead per spawn

## Blocked On

Track upstream progress:
- [Issue #2604](https://github.com/openai/codex/issues/2604) - 261+ upvotes, main tracking issue
- [Issue #8664](https://github.com/openai/codex/issues/8664) - spawn_subagents/chain_subagents proposal
- [PR #3655](https://github.com/openai/codex/pull/3655) - multi subagent orchestration (in progress)

## When Unblocked

### Tasks
1. Update `scripts/install-codex.sh` to install `agents/` directory
2. Create agent format converter (Claude Code → Codex)
3. Update flow-next skills to detect Codex + use spawn_agent
4. Test parallel scout execution during `/flow-next:plan`
5. Document Codex-specific limitations

### Scope
- Planning phase scouts (practice, docs, repo, github)
- Worker agent for `/flow-next:work`
- Skip: memory-scout, plan-sync (require Claude Code features)

## Quick commands
- `./scripts/install-codex.sh flow-next` - current install
- `.flow/bin/flowctl --help` - CLI reference

## Acceptance
- [ ] Custom agent definitions work in Codex
- [ ] Parallel scouts during planning phase
- [ ] Worker agent for task execution
- [ ] Install script copies agents/
- [ ] Documented in README

## References
- Codex changelog: https://developers.openai.com/codex/changelog/
- Codex skills: https://developers.openai.com/codex/skills/
- Agents SDK: https://developers.openai.com/codex/guides/agents-sdk/
- Orchestrator template: `codex-rs/core/templates/agents/orchestrator.md`
- AgentRole enum: `codex-rs/core/src/agent/role.rs`
- codex-subagents-mcp: https://github.com/leonardsellem/codex-subagents-mcp
