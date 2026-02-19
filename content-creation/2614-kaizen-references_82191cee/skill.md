# Kaizen Plugin — Reference Materials

**Date fetched:** 2026-02-17

## User-Provided Analysis Dimensions (beyond original scope)

1. **Shortest path to outcome** — where is energy and steps being wasted?
2. **Red herrings** — distractions that get chased often but are always just a distraction from the actual solution
3. **System process interruptions** — what system processes interrupted correct paths?
4. **Missing hooks** — what hooks could mould behavior and remove AI overhead by providing details automatically and automating correction?
5. **MCP-based DuckDB for transcript querying** — MotherDuck MCP server provides SQL analytics directly to AI assistants

## Claude Code Hooks — Key Facts for Kaizen Analysis

Source: <https://code.claude.com/docs/en/hooks-guide.md>, <https://code.claude.com/docs/en/hooks.md> (accessed 2026-02-17)

### Hook Events (Full List)

- `SessionStart` — session begins or resumes (matchers: startup, resume, clear, compact)
- `UserPromptSubmit` — user submits prompt, before processing
- `PreToolUse` — before tool call executes (can block, allow, deny, modify input)
- `PermissionRequest` — permission dialog appears
- `PostToolUse` — after tool call succeeds
- `PostToolUseFailure` — after tool call fails
- `Notification` — Claude sends notification (matchers: permission_prompt, idle_prompt, auth_success, elicitation_dialog)
- `SubagentStart` — subagent spawned
- `SubagentStop` — subagent finishes
- `Stop` — Claude finishes responding
- `TeammateIdle` — agent team teammate about to go idle
- `TaskCompleted` — task being marked completed
- `PreCompact` — before context compaction (matchers: manual, auto)
- `SessionEnd` — session terminates

### Hook Types

- `command` — shell command (receives JSON stdin, returns via exit code + stdout/stderr)
- `prompt` — single-turn LLM evaluation (Haiku default, returns `{ok: true/false, reason}`)
- `agent` — multi-turn subagent with tool access (up to 50 turns, returns same format)

### Hook Input (Common Fields)

```json
{
  "session_id": "abc123",
  "transcript_path": "/path/to/transcript.jsonl",
  "cwd": "/path/to/project",
  "permission_mode": "default",
  "hook_event_name": "PreToolUse",
  "tool_name": "Bash",
  "tool_input": { "command": "..." }
}
```

### PreToolUse Decision Control

```json
{
  "hookSpecificOutput": {
    "hookEventName": "PreToolUse",
    "permissionDecision": "deny",
    "permissionDecisionReason": "Use Glob instead of Bash ls",
    "updatedInput": { "command": "modified command" },
    "additionalContext": "extra info for Claude"
  }
}
```

Decisions: `allow` (bypass permission), `deny` (block + reason to Claude), `ask` (prompt user)

### Key Capabilities for Kaizen

1. **PreToolUse hooks can deny and redirect** — when analysis identifies anti-patterns (e.g., Bash ls), a hook can deny and tell Claude to use the correct tool
2. **PostToolUseFailure hooks** — can inject corrective context after failures
3. **SubagentStart hooks** — can inject context into subagents (e.g., "write research to files, not messages")
4. **SubagentStop hooks** — can validate subagent output quality before accepting
5. **Stop hooks** — can verify completeness before allowing Claude to stop
6. **SessionStart (compact matcher)** — can re-inject critical context lost during compaction
7. **TeammateIdle hooks** — can enforce quality gates before teammates go idle
8. **TaskCompleted hooks** — can enforce completion criteria
9. **PreCompact hooks** — can save state before compaction

### Hooks in Skills and Agents (Frontmatter)

Hooks can be defined in skill/agent YAML frontmatter — scoped to component lifecycle:

```yaml
---
name: secure-operations
hooks:
  PreToolUse:
    - matcher: "Bash"
      hooks:
        - type: command
          command: "./scripts/security-check.sh"
---
```

## Headless Mode / Agent SDK

Source: <https://code.claude.com/docs/en/headless.md> (accessed 2026-02-17)

- `claude -p "prompt"` runs non-interactively
- `--output-format json` returns structured JSON with session ID
- `--output-format stream-json` for real-time streaming
- `--allowedTools` for auto-approving specific tools
- `--json-schema` for structured output conforming to schema
- Can continue conversations with `--continue` or `--resume <session_id>`
- `--append-system-prompt` adds instructions while keeping defaults

### Relevance to Kaizen

Headless mode enables automated transcript analysis pipelines:

```bash
# Analyze a transcript with Claude itself
claude -p "Analyze this session for anti-patterns" \
  --allowedTools "Read,Grep,Glob" \
  --output-format json \
  --json-schema '{"type":"object","properties":{"anti_patterns":{"type":"array"}}}'
```

## MotherDuck MCP Server for DuckDB

Source: <https://github.com/motherduckdb/mcp-server-motherduck> (accessed 2026-02-17)

### What It Is

MCP server that gives AI assistants direct SQL access to DuckDB databases. Supports local files, in-memory, S3, and MotherDuck cloud.

### Tools Provided

- `execute_query` — execute SQL (DuckDB dialect)
- `list_databases` — list all databases
- `list_tables` — list tables and views
- `list_columns` — list columns of a table/view
- `switch_database_connection` — switch to different database (requires --allow-switch-databases)

### Claude Code Integration

```bash
claude mcp add --scope user duckdb --transport stdio -- \
  uvx mcp-server-motherduck --db-path :memory: --read-write --allow-switch-databases
```

### Relevance to Kaizen

Instead of writing custom Python scripts for DuckDB queries, the kaizen plugin could include an MCP server configuration that lets Claude query transcript data directly via SQL during analysis:

```sql
SELECT
  json_extract_string(line, '$.type') as event_type,
  json_extract_string(line, '$.message.content[0].name') as tool_name,
  COUNT(*) as frequency
FROM read_ndjson_auto('/path/to/transcripts/*.jsonl')
WHERE json_extract_string(line, '$.type') = 'assistant'
GROUP BY event_type, tool_name
ORDER BY frequency DESC;
```

This gives the analysis agent SQL querying as a native tool rather than requiring custom scripts.
