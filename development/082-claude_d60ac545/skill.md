# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Claudest is a curated Claude Code plugin marketplace containing five plugins: **claude-memory** (conversation memory with full-text search and context injection), **claude-utilities** (convert-to-markdown via ezycopy), **claude-skills** (skill authoring and repair), **claude-coding** (git workflows and CLAUDE.md maintenance), and **claude-thinking** (structured thinking tools). There is no build system or package manager — plugin runtime is Python 3.7+ stdlib-only. Tests use pytest with hypothesis (dev dependencies only).

## Development Commands

```bash
# Re-import all conversations into the memory DB
python3 plugins/claude-memory/hooks/import_conversations.py

# Import with stats output
python3 plugins/claude-memory/hooks/import_conversations.py --stats

# Test context injection manually
echo '{"source":"startup","session_id":"test","cwd":"/some/path"}' | python3 plugins/claude-memory/hooks/memory-context.py

# Test session sync manually
echo '{"session_id":"<uuid>"}' | python3 plugins/claude-memory/hooks/sync_current.py
```

No `--force` flag exists for import — delete `~/.claude-memory/conversations.db` and re-run to reimport from scratch.

## Architecture

### Plugin Structure

Each plugin follows the Claude Code plugin convention: `.claude-plugin/plugin.json` for metadata, `hooks/` for lifecycle hooks, `commands/` for slash commands, and `skills/` for auto-triggered capabilities.

### claude-memory Hook Lifecycle

On **SessionStart**, two hooks fire in order: `memory-setup.py` creates the `~/.claude-memory/` directory and kicks off a background import if the DB doesn't exist, then `memory-context.py` (on `startup|clear` events only) queries recent sessions and injects them as additional context via `hookSpecificOutput`. On **Stop**, `memory-sync.py` writes hook input to a temp file and spawns `sync_current.py --input-file` in the background to incrementally sync the current session to the DB without blocking shutdown. All hooks are Python (no bash) for cross-platform compatibility.

### Database (v3 Schema)

SQLite at `~/.claude-memory/conversations.db` with WAL mode and 5s busy_timeout for concurrent access safety. The key design: messages are stored once per session (deduped), and branches (from conversation rewinds) are tracked via a many-to-many `branch_messages` table. Full-text search uses a cascade: FTS5 (best, BM25 ranking) → FTS4 (MATCH + snippet, no BM25) → LIKE fallback (no stemming). Schema is split into `SCHEMA_CORE` (tables/indexes) and `SCHEMA_FTS5`/`SCHEMA_FTS4` variants, applied conditionally based on `detect_fts_support()`. Auto-migrated on connection — if the schema version is outdated, the DB is deleted and recreated.

Tables: `projects`, `sessions`, `branches`, `messages`, `branch_messages`, `import_log`, `messages_fts` (virtual), `branches_fts` (virtual).

### Shared Code

`plugins/claude-memory/skills/recall-conversations/scripts/memory_lib/` is the shared utility package used by all hooks and skill scripts. It is split into four focused modules: `db.py` (database connection, schema, settings, logging), `content.py` (message content extraction and tool detection), `parsing.py` (JSONL parsing, branch detection via UUID parent chain analysis, metadata extraction), and `formatting.py` (session formatting, time/path utilities). The `extract_text_content()` function in `content.py` returns a 4-tuple `(text, has_tool_use, has_thinking, tool_summary_json)` — tool markers are never materialized into stored text; instead tool counts are stored as compact JSON in `messages.tool_summary`.

### Session Selection Algorithm

`memory-context.py:select_sessions` iterates recent sessions (excluding current session and subagents), skips those with `exchange_count <= 1`. For 2-exchange sessions it appends and keeps looking (up to `max_context_sessions`). For sessions with >2 exchanges it appends and stops. This means multiple sessions are only injected when the most recent ones are very short.

## Conventions

Commit messages use conventional commits: `feat(memory):`, `fix(memory):`, `chore(memory):`, `docs:`, `refactor(memory):`. Version is tracked in two places that must stay in sync: each plugin's `.claude-plugin/plugin.json` and the root `.claude-plugin/marketplace.json`. Always bump version before pushing changes.

Settings are hardcoded defaults in `memory_lib/db.py:DEFAULT_SETTINGS` (the YAML settings file was removed since PyYAML is not stdlib and settings were silently ignored for most users).

Skill descriptions in SKILL.md frontmatter should be short and focused — verbose descriptions pollute the agent's context window since they're loaded whenever the skill triggers.
