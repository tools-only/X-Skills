# Design: hooks-guide skill for plugin-creator

**Date:** 2026-02-27
**Status:** Approved
**Plugin:** `plugins/plugin-creator`

---

## Problem

The existing hooks skills in `plugin-creator` cover Claude Code hooks only (Node.js `.cjs`,
plugin-level, settings-level). The AI assistant ecosystem has expanded — GitHub Copilot, Cursor,
Windsurf, Amp, and others now publish hook systems. These docs are human-facing and change
frequently. There is no single skill that:

- Covers hooks across multiple platforms with verified, AI-optimised content
- Documents hook authoring in both Node.js CJS and Python
- Includes inline-agent hooks/MCP (new Claude Code sub-agents convention)
- Stays current via a fetch-and-transform pipeline

---

## Decision: Option A

Single `hooks-guide` skill with a navigation SKILL.md and per-topic/per-platform reference files
in `references/`. A fetch script pulls raw human-facing docs per platform and runs them through
the `rwr:doc-to-skill` pipeline to produce AI-facing reference content.

---

## Skill location

```
plugins/plugin-creator/skills/hooks-guide/
├── SKILL.md                          ← navigation hub, decision flowcharts
├── scripts/
│   └── fetch-and-transform-hooks-docs.sh   ← fetch + rwr:doc-to-skill pipeline
└── references/
    ├── common-schema.md              ← shared concepts: events, matchers, exit codes, JSON I/O
    ├── claude-code.md                ← Claude Code hooks (sourced from fetched hooks-doc.md)
    ├── inline-agent-hooks.md         ← inline hooks/MCP in agent frontmatter (sub-agents doc)
    ├── github-copilot.md             ← GitHub Copilot coding agent hooks
    ├── hooks-cjs.md                  ← Node.js CommonJS authoring guide
    ├── hooks-python.md               ← Python authoring guide
    ├── best-practices.md             ← cross-platform conventions and anti-patterns
    └── platform-coverage.md          ← known platforms, fetch URLs, coverage status
```

---

## SKILL.md design

**Frontmatter:**

```yaml
---
description: >
  Cross-platform hooks reference for AI coding assistants — Claude Code, GitHub Copilot,
  Cursor, Windsurf, Amp. Covers hook authoring in Node.js CJS and Python, per-platform
  event schemas, inline-agent hooks and MCP in agent frontmatter, common JSON I/O,
  exit codes, best practices, and fetch scripts to refresh docs from official sources.
  Use when writing, reviewing, or debugging hooks for any AI assistant.
allowed-tools: Read, Grep, Glob, Bash, Write, Edit
---
```

**Body structure:**

1. Quick-start decision flowchart — which platform? which language?
2. Common concepts section (events, matchers, stdin JSON, exit codes) — load `common-schema.md`
3. Per-platform routing — points to the correct reference file
4. Language guide routing — CJS vs Python
5. Inline-agent hooks routing — `inline-agent-hooks.md`
6. Refresh instructions — how to run the fetch script

---

## Reference file content sources

| File | Primary source | Transform method |
|------|---------------|-----------------|
| `claude-code.md` | `.claude/plan/plugin-creator-hooks/hooks-doc.md` (already fetched) | `rwr:doc-to-skill` |
| `inline-agent-hooks.md` | `.claude/plan/plugin-creator-hooks/sub-agents-doc.md` (already fetched) | `rwr:doc-to-skill` — hooks + mcpServers + skills sections only |
| `github-copilot.md` | `.claude/plan/plugin-creator-hooks/github-copilot-hooks-doc.md` (already fetched) | `rwr:doc-to-skill` |
| `hooks-cjs.md` | Derived from existing `hooks-core-reference` + `hooks-patterns` | Authored — no external fetch needed |
| `hooks-python.md` | Authored from Claude Code hooks doc stdin/exit-code spec | Authored — Python stdlib examples |
| `common-schema.md` | Synthesised from all platform docs | Authored — cross-platform normalisation |
| `best-practices.md` | Synthesised from all platform docs | Authored |
| `platform-coverage.md` | Maintained manually | Static — lists platforms, URLs, last-verified dates |

---

## Fetch-and-transform pipeline

The script `scripts/fetch-and-transform-hooks-docs.sh` does the following for each platform:

1. Fetch raw human-facing doc via `curl` to a temp file in `.claude/plan/plugin-creator-hooks/`
2. Invoke `rwr:doc-to-skill` (via `claude -p`) on the temp file, writing output to
   `skills/hooks-guide/references/<platform>.md`
3. Log fetch timestamp and source URL into `platform-coverage.md`
4. Exit non-zero if any platform fetch fails, but continue other platforms (graceful partial)

**Platforms and URLs for initial implementation:**

| Platform | URL | Status |
|----------|-----|--------|
| Claude Code | `https://code.claude.com/docs/en/hooks.md` | Verified — fetched |
| Claude Code sub-agents | `https://code.claude.com/docs/en/sub-agents.md` | Verified — fetched |
| GitHub Copilot | `https://docs.github.com/en/copilot/how-tos/use-copilot-agents/coding-agent/use-hooks.md` | Verified — fetched |
| Cursor | `https://docs.cursor.com/context/rules` | To verify at fetch time |
| Windsurf | `https://docs.windsurf.com/windsurf/memories` | To verify at fetch time |
| Amp | `https://ampcode.com/docs` | To verify at fetch time |

Platforms with no hook concept return 404 or no hooks section — script logs this to
`platform-coverage.md` as "no hooks documented" and skips reference file update.

---

## Changes to existing files

### `agents/hook-creator.md`

1. Add `plugin-creator:hooks-guide` to the `skills:` preload list in frontmatter
2. Add new branch to Scope Determination flowchart:
   - "Scoped to one agent only" → "Inline agent frontmatter hooks" → `inline-agent-hooks.md`
3. Add five new events to Event Selection flowchart:
   - `TeammateIdle`, `TaskCompleted`, `ConfigChange`, `WorktreeCreate`, `WorktreeRemove`
4. Add `SubagentStart` `additionalContext` injection pattern to Phase 2 script templates
5. Update Sources section with new doc fetch dates

### `skills/claude-hooks-reference-2026/SKILL.md`

Add `plugin-creator:hooks-guide` to the skills loaded by the umbrella, so when the umbrella
is activated the broader cross-platform guide is also available.

---

## What is NOT changing

- `skills/hooks-core-reference/`, `hooks-io-api/`, `hooks-patterns/` — unchanged, still the
  Claude Code deep-reference trio
- `skills/claude-hooks-reference-2026/` — unchanged structure, just gains one more skill in its
  load list
- No new agent — `hook-creator` is updated, not replaced or duplicated
- No Cursor/Windsurf content fabricated — fetch script attempts URLs, logs result, skips if
  no hooks content found

---

## Validation gates

After implementation:

```bash
uv run plugins/plugin-creator/scripts/plugin_validator.py plugins/plugin-creator
```

All new and modified files must pass. Token budget per reference file: stay under
`TOKEN_WARNING_THRESHOLD` defined in `plugin_validator.py`.

---

## Out of scope

- Hook testing infrastructure (separate backlog item)
- Cursor MDC / rules system (different concept from hooks — separate skill if needed)
- Windsurf memories system (different concept — not hooks)
- PowerShell hook scripts (Copilot supports it; out of scope for v1)
