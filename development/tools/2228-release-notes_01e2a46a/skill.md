# Release Notes

## v2.1 — Auto Mode + Deterministic Enforcement + Portability

**Date**: 2026-02-12

### Why This Release

Two problems drove this release:

1. **Hook enforcement was too aggressive for non-autonomous users.** The stop-validator checkpoint gate activated for *every* session that made code changes — even casual interactive use. Users who didn't want autonomous workflows were blocked from stopping until they wrote a full checkpoint with `key_insight`, `search_terms`, `linters_pass`, etc. The fix: make the entire enforcement system opt-in via `/auto`.

2. **The toolkit was not portable.** Hooks contained hardcoded project paths, doc counts, and file patterns from the original development repo. Moving the toolkit to any other project caused silent degradation — hooks that checked for specific files would silently skip, scoring functions would return zero for every event.

### What Changed

#### Auto Mode Toggle (`/auto`)

The entire hook enforcement system is now opt-in. Run `/auto` (or `/auto on`) to activate autonomous mode; `/auto off` to deactivate.

| Hook Behavior | Auto OFF (default) | Auto ON |
|---------------|-------------------|---------|
| Stop-validator checkpoint | Not required | Required for code-changing sessions |
| Auto-approve tool calls | Disabled | Enabled |
| Deploy enforcer | Disabled | Active |
| State file guard | Disabled | Active |
| Skill continuation | Disabled | Active |
| Memory capture on stop | Opportunistic (if checkpoint exists) | Active (always) |

**Key change**: `stop-validator.py:requires_checkpoint()` now returns `False` when autonomous mode is not active. Previously it returned `True` for any session with code changes, blocking all users.

Skills that activate autonomous mode (`/melt`, `/repair`, `/forge`, etc.) still work exactly as before — they create `autonomous-state.json` which activates the hooks. The difference is that *without* a skill or `/auto`, hooks are completely passive.

```bash
# Toggle on
/auto           # or /auto on

# Toggle off
/auto off

# Check status
/auto status
```

#### Deterministic Quality Enforcement (`/enforce`)

New skill that researches the repo's tech stack, creates strict linter configs, and fixes every violation across the entire codebase.

- Ships strict reference configs for ruff and ESLint via `config/references/` (used when the project has no config)
- Enforces structural limits: files < 400 LOC, functions < 80 LOC
- Integrated into stop-validator: `_linter.py` module runs deterministic linters as a checkpoint gate
- System-detected code changes (`_code_changes_detected`) replace self-reported `code_changes_made` boolean — the model cannot bypass linter gates by lying

#### Stop-Validator Hardening

The stop-validator received its most significant upgrade since creation — 22 new tests, anti-fabrication gates, and workflow phase enforcement:

| Gate | What it prevents |
|------|-----------------|
| Test execution proof | Cross-references `tool-usage-log.json` for Bash usage when checkpoint claims command_output test results — blocks fabricated test results |
| Workflow phase proof | When `autonomous-state.json` has a `workflow` field, verifies browser discovery tools were used *before* code editing tools |
| Phased work detection | Catches agents writing `is_job_complete: true` after completing one phase of multi-phase work |
| Acceptance criteria enforcement | Tamper-resistant — `acceptance-criteria.json` is written by hooks only, protected by state-file-guard |
| Inherited dirty state | Pre-existing uncommitted changes no longer trigger false blocking — uses dirty file list instead of diff hash |
| Terminal noise reduction | Blocking messages compacted to essential information only |

#### Portability

All project-specific references removed. The toolkit now works in any repository without modification:

- Stripped hardcoded file paths, doc counts, and project names from all hooks
- Dynamic discovery via filesystem globs replaces hardcoded maps
- Portable quality configs shipped in `config/references/` (not generated per-repo)
- 6 portability leaks found and fixed by heavy verification agents
- `/audit` command added for installation health checks

#### New Skills

| Skill | Purpose |
|-------|---------|
| `/auto` | Toggle hook enforcement on/off |
| `/forge` | Experience-first autonomous debugging (discover via browser, then fix) |
| `/enforce` | Deterministic quality enforcement across entire codebase |
| `/loom` | Browser walkthrough recording with AI voiceover |
| `/plugin` | SELD relay plugin authoring lifecycle |
| `/essay` | Publication-quality essay writing with 5-agent analysis |

#### Video Skill Rewrite

Complete rewrite to Feb 2026 SOTA models and a directing paradigm:

- Default engine: Kling 3.0 V3 Pro (cinematic quality) with character locking via Elements
- Multi-shot generation is default — one prompt produces connected scenes
- Claude acts as Creative Director, not prompt engineer
- Removed Remotion dependency — assembly via DaVinci Resolve / CapCut / FFmpeg

#### Memory System Improvements

- Inverted entity index built at write time (amortized O(1) per entity) — candidate pool goes from 25% to 46% corpus reach
- 2-signal scoring (entity 60% + recency 40%) replaces 3-signal (entity + recency + quality) — quality signal was saturating at 1.0
- Synonym expansion in index query — semantically equivalent terms (`race-condition` / `concurrency`) now match
- Implicit memory citation detection — scans `what_was_done` and `key_insight` for `m1`-`m9` patterns

#### Modular CLAUDE.md

`CLAUDE.md` split into `~/.claude/rules/toolkit-*.md` files:

| File | Content |
|------|---------|
| `toolkit-skills.md` | Skill routing table and parallelization strategy |
| `toolkit-git.md` | Git operations in autonomous mode |
| `toolkit-search.md` | Exa MCP and QMD search preferences |

This eliminates merge conflicts when the toolkit updates — rules files are managed automatically, user preferences stay in `CLAUDE.md`.

#### Debugging Infrastructure

Three tiers of debugging infrastructure added:

- **Tier 1**: Hook health validation at session start, verification artifact cross-check, utility feedback loop
- **Tier 2**: Cross-project memory sharing, crash recovery checkpoints, debugging-aware scoring
- **Tier 3**: Bash error advisor (pattern matching for common failures), autonomous health monitor, deploy rollback

### Configuration Changes

- `CLAUDE_CODE_MAX_OUTPUT_TOKENS` reduced from 128K to 64K (cost optimization, no quality impact)
- `autonomous-state.json` removed from state-file-guard's `PROTECTED_FILES` (skills and `/auto` need to write it)
- Agent-type Stop hook removed from settings.json (deterministic Python validation is sufficient)
- `_linter.py` added to shared module health checks in session-init

### New Files

| File | Purpose |
|------|---------|
| `config/commands/auto.md` | `/auto` toggle command definition |
| `config/hooks/_linter.py` | Shared linter module for deterministic quality gates |
| `config/references/ruff-strict.toml` | Strict ruff config used when project has none |
| `config/references/eslint-strict.config.mjs` | Strict ESLint config used when project has none |
| `config/rules/toolkit-skills.md` | Skill routing (split from CLAUDE.md) |
| `config/rules/toolkit-git.md` | Git conventions (split from CLAUDE.md) |
| `config/rules/toolkit-search.md` | Search tool preferences (split from CLAUDE.md) |

### Migration Guide

**From v2.0**:

1. Re-run the installer: `./scripts/install.sh`
2. Restart Claude Code
3. Hooks are now passive by default. To get v2.0 behavior, run `/auto on` at the start of autonomous sessions
4. Skills (`/melt`, `/repair`, `/forge`, etc.) still activate autonomous mode automatically — no change needed for skill-driven workflows

**For teams with non-autonomous users**: No action needed. The default experience is now non-blocking. Users who want autonomous enforcement opt in with `/auto`.

### Known Limitations

- `/auto status` reads `autonomous-state.json` directly — if the file was created by a skill, the `mode` field will show the skill name (e.g., `"melt"`) rather than `"auto"`
- Deterministic linter gates require the linter binary to be installed (ruff, eslint, etc.) — if missing, the gate is skipped with a warning
- State-file-guard no longer protects `autonomous-state.json` — an agent could theoretically delete it to bypass enforcement (mitigated by the hook still checking for the file's existence)

---

## v2.0 — Opus 4.6 + Agent Swarms + Complementary Memory

**Date**: 2026-02-07

### Why This Release

Claude Code has evolved rapidly. Two major platform changes drove this release:

1. **Opus 4.6** is significantly more capable than earlier models — it follows instructions more naturally, reasons better about multi-step tasks, and doesn't need the rigid guardrails (mandatory agent counts, ASCII checkpoint diagrams, forced verification steps) that earlier models required. We've refactored all 26 skills from prescriptive procedures into natural capability modes.

2. **Agent Teams** (experimental) enables true parallel agent swarms — multiple Claude instances coordinating via shared task lists, direct messaging, and team-scoped context. This replaces the previous pattern of sequential `Task()` calls with genuine concurrent collaboration.

3. **Native MEMORY.md** (Claude Code v2.1.32+) adds built-in project memory that loads unconditionally every session. Our custom compound memory system needed to integrate with it rather than compete for context budget.

### What Changed

#### Agent Swarms (Agent Teams)

Skills now leverage Agent Teams for parallel work:

- **`/melt`** spawns planning swarms (First Principles + AGI-Pilled + domain experts), then implementation agents for independent work items
- **`/heavy`** spawns 3-5 analysis agents for simultaneous multi-perspective debate
- **`/burndown`** spawns detection agents that scan different codebase aspects concurrently

Teams coordinate through shared task lists (`~/.claude/tasks/`) with blocking dependencies, auto-assignment, and structured messaging. When a teammate completes a task, it's automatically communicated to the team lead.

**How to enable**:

Add this env var to your `~/.claude/settings.json` (global) or `config/settings.json` (repo):

```json
{
  "env": {
    "CLAUDE_CODE_EXPERIMENTAL_AGENT_TEAMS": "1"
  }
}
```

The install script sets this automatically. Restart Claude Code after changing settings.

> **Note**: Agent Teams is an experimental Claude Code feature. It's enabled via environment variable, not a CLI flag. Requires Claude Code v2.1.32+.

**Without Agent Teams**: All skills still work — they fall back to sequential `Task()` subagent calls. You lose shared task lists and team coordination, but the core fix-verify loop is identical.

#### Complementary Memory (Native + Custom)

The compound-context-loader now auto-detects Claude Code's native MEMORY.md and adjusts its behavior:

| Feature | Before | After |
|---------|--------|-------|
| Compound memory budget | Fixed 8000 chars | Dynamic: 4500 with MEMORY.md, 8000 without |
| Native memory awareness | None | Auto-detected, dedup guard active |
| Content deduplication | None cross-system | Events >60% overlapping MEMORY.md are skipped |
| Memory promotion | Manual only | Utility-based script promotes top events |

**How it works**:

1. At session start, the hook checks for `~/.claude/projects/{encoded-path}/memory/MEMORY.md`
2. If found, it reduces compound memory budget from 8K to 4.5K chars (native memory uses ~4-6K)
3. A dedup guard skips compound events whose LESSON content is already well-documented in MEMORY.md
4. High-utility events (cited/injected ratio >= 30%) can be promoted to MEMORY.md

**Promotion workflow**:

```bash
# See which events qualify
python3 config/scripts/promote-to-memory-md.py --dry-run

# Promote top 3 to MEMORY.md
python3 config/scripts/promote-to-memory-md.py
```

Promoted events are tracked in a sidecar file (`promoted-events.json`) so they aren't re-promoted.

**Zero regression risk**: If MEMORY.md detection fails, the system falls back to the 8000-char standalone budget — identical to previous behavior.

#### Opus 4.6 Skill Optimization

All 26 skills were refactored from prescriptive procedures to natural capability modes:

- Removed mandatory agent counts (Opus 4.6 naturally spawns the right number)
- Removed forced ASCII checkpoint diagrams with incorrect fields
- Removed hook enforcement labels that confused more than helped
- Replaced with natural encouragement that the model follows better
- Single universal checkpoint schema from `stop-validator.py` is the only truth

**Key insight**: Skills are MODES that enable capabilities, not PROCEDURES that prescribe steps. When the model is capable enough, removing guardrails produces cleaner execution.

#### Hook System Consolidation

Unified from 24 hooks to 13, reducing ~4500 lines to ~1500:

- Merged duplicate event handlers
- Removed legacy status-file hooks (replaced by Claude Code's native Tasks)
- Combined overlapping functionality into shared modules (`_memory.py`, `_common.py`, `_session.py`)

### Configuration Changes

#### settings.json

```json
{
  "env": {
    "CLAUDE_CODE_MAX_OUTPUT_TOKENS": "128000",
    "ENABLE_TOOL_SEARCH": "auto",
    "CLAUDE_CODE_EXPERIMENTAL_AGENT_TEAMS": "1"
  },
  "permissions": {
    "defaultMode": "dontAsk"
  },
  "alwaysThinkingEnabled": true
}
```

**New env vars**:
- `CLAUDE_CODE_EXPERIMENTAL_AGENT_TEAMS=1` — enables Agent Teams (swarm mode)
- `ENABLE_TOOL_SEARCH=auto` — defers MCP tool loading until needed (saves 85-95% context tokens)

### New Files

| File | Purpose |
|------|---------|
| `config/scripts/promote-to-memory-md.py` | Promote high-utility compound memories to native MEMORY.md |
| `~/.claude/projects/.../memory/MEMORY.md` | Native project memory seed (created by toolkit) |

### Migration Guide

**From v1.x (pre-Opus 4.6)**:

1. Re-run the installer: `./scripts/install.sh`
2. Restart Claude Code
3. Agent Teams is auto-enabled. To disable: remove `CLAUDE_CODE_EXPERIMENTAL_AGENT_TEAMS` from settings.json
4. MEMORY.md is auto-detected. No action needed — the toolkit seeds it if empty

**Manual Agent Teams setup** (if not using install script):

```bash
# Add to ~/.claude/settings.json
python3 -c "
import json
from pathlib import Path
p = Path.home() / '.claude' / 'settings.json'
s = json.loads(p.read_text()) if p.exists() else {}
s.setdefault('env', {})['CLAUDE_CODE_EXPERIMENTAL_AGENT_TEAMS'] = '1'
p.write_text(json.dumps(s, indent=2))
print('Agent Teams enabled. Restart Claude Code.')
"
```

### Known Limitations

- Agent Teams is experimental — teammates occasionally go idle between turns (this is normal; they wake on message)
- For single-agent read-only analysis, a direct Explore agent is more efficient than a full team
- MEMORY.md path encoding assumes Claude Code's convention (replace `/` and `_` with `-`); a glob fallback handles edge cases
