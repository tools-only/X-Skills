---
name: session-historian
description: 'Look up prior Claude Code sessions when context is lost or forgotten. Use when asked "what did we do before?", "what happened in the last session?", "I forgot what we were working on", "find what I told you about X", or any request to recall past conversation history, prior decisions, experiments, or outcomes. Searches raw JSONL transcripts from ~/.claude/projects/ via DuckDB index. Returns verbatim user messages and summarizes AI actions and sub-agent outcomes. Summaries cached at ~/.claude/kaizen/session-summaries/.'
---

# Session Historian

Look up what happened in prior Claude Code sessions. Index and search raw JSONL transcripts. Cache session summaries.

## Script

```text
.claude/skills/session-historian/scripts/session_query.py
```

Run directly (shebang + PEP 723, no manual `uv run` prefix needed):

```bash
.claude/skills/session-historian/scripts/session_query.py list
.claude/skills/session-historian/scripts/session_query.py messages <session-id>
.claude/skills/session-historian/scripts/session_query.py search "experiment worktree"
.claude/skills/session-historian/scripts/session_query.py show <session-id>
.claude/skills/session-historian/scripts/session_query.py index
```

## Workflow: "I forgot what happened"

1. **List recent sessions** to find the relevant one:

   ```bash
   .claude/skills/session-historian/scripts/session_query.py list --limit 10
   ```

2. **Read verbatim user messages** (fastest way to reconstruct context):

   ```bash
   .claude/skills/session-historian/scripts/session_query.py messages <session-id-prefix>
   ```

3. **Search for specific topics** across all raw files:

   ```bash
   .claude/skills/session-historian/scripts/session_query.py search "experiment" --limit 10
   ```

4. **Generate and cache a session summary** (when detail is needed):
   - Run `show <session-id>` to get session metadata and file path
   - Read the raw JSONL file
   - Produce a structured summary following the template below
   - Write to `~/.claude/kaizen/session-summaries/<session-id>.md`
   - Run `mark-summarized <session-id>` to flag in the index

## Summary Template

When generating a session summary, write to `~/.claude/kaizen/session-summaries/<session-id>.md`:

```markdown
---
source_type: claude-session
source_path: <file_path>
session_id: <session_id>
project: <project_name>
date_range: <first_ts> → <last_ts>
user_message_count: <N>
assistant_turns: <N>
confidence: high|medium|low
generated_at: <ISO timestamp>
---

## Summary

[BLUF: 2-3 sentences on what this session accomplished]

## User Messages (verbatim)

[Numbered list of every user message, exactly as written, with timestamp]

1. [2026-02-21T15:40] "First ensure that .claude/worktrees is gitignored..."
2. [2026-02-21T16:51] "Ah. Do you not remember the experiment..."

## Actions Taken

[Bullet list of what the AI did — file edits, commands run, commits made]

- Created .claude/worktrees directory and added to .gitignore
- Spawned sub-agent a0263e5 to read workflow docs → produced project_workflow.draft.md
- ...

## Sub-Agent Tasks and Outcomes

| Agent ID | Task | Outcome |
|----------|------|---------|
| a0263e5  | Read docs, produce workflow draft | Produced project_workflow.draft.md |
| a5d97a8  | Execute work-backlog-item workflow | Partial — blocked at AskUserQuestion node |

## Key Findings / Decisions

[Important decisions made, things discovered, conclusions reached]

## What Was NOT Completed

[Work started but not finished, blocked items]

## Sources

File: <file_path>
Indexed: <DB_PATH>
```

## Fidelity Rules (from summarizer skill)

- **Read before summarizing**: Read the actual JSONL file, do not guess from metadata.
- **Verbatim user messages**: Copy exactly — never paraphrase.
- **Preserve counts**: "3 sub-agents spawned, 2 succeeded, 1 blocked" not "most sub-agents succeeded."
- **Distinguish absence**: "Not mentioned in transcript" not "didn't happen."
- **Search raw files**: `search` command reads raw JSONL, not summaries. Never search summaries for facts.

## Index Location

```text
~/.claude/kaizen/session-index.duckdb    ← DuckDB index (sessions + user_messages tables)
~/.claude/kaizen/session-summaries/      ← Cached AI-generated summaries
```

The index is built on first `list` or explicitly with `index`. It is NOT automatically kept current — re-run `index` or pass `--rebuild` to pick up new sessions.

## JSONL Schema Reference

Records in `~/.claude/projects/<project-slug>/<session-id>.jsonl`:

```json
{"type": "user", "sessionId": "uuid", "timestamp": "ISO8601",
 "message": {"content": "string or [{type,text}]"}}
{"type": "assistant", "sessionId": "uuid", "timestamp": "ISO8601",
 "message": {"content": [{...}]}}
```

Records with `toolUseResult` key are tool results — skip for user message extraction.
Tool use blocks have `type: "tool_use"` inside `message.content` array.
