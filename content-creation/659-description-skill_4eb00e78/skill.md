---
description: Hook JSON input/output API reference — what data hooks receive via stdin and what JSON they can return to control Claude Code behavior. Use when writing hook scripts, checking exit code behavior, building JSON output for PreToolUse permissions, or understanding event-specific input schemas.
user-invocable: true
---

# Claude Code Hooks — I/O API Reference (January 2026)

JSON schemas for hook stdin input and stdout output per event. For hook system fundamentals, activate `Skill(command: "plugin-creator:hooks-core-reference")`. For working examples, activate `Skill(command: "plugin-creator:hooks-patterns")`.

---

## Hook Input (JSON via stdin)

### Common Fields (All Events)

```json
{
  "session_id": "abc123",
  "transcript_path": "/path/to/session.jsonl",
  "cwd": "/current/working/directory",
  "permission_mode": "default",
  "hook_event_name": "PreToolUse"
}
```

**`permission_mode` values**: `"default"`, `"plan"`, `"acceptEdits"`, `"dontAsk"`, `"bypassPermissions"`

### PreToolUse Input

```json
{
  "session_id": "abc123",
  "transcript_path": "/path/to/session.jsonl",
  "cwd": "/path/to/project",
  "permission_mode": "default",
  "hook_event_name": "PreToolUse",
  "tool_name": "Bash",
  "tool_input": {
    "command": "psql -c 'SELECT * FROM users'",
    "description": "Query the users table",
    "timeout": 120000
  },
  "tool_use_id": "toolu_01ABC123"
}
```

**Tool-specific `tool_input` fields:**

| Tool    | Fields                                                   |
| ------- | -------------------------------------------------------- |
| `Bash`  | `command`, `description`, `timeout`, `run_in_background` |
| `Write` | `file_path`, `content`                                   |
| `Edit`  | `file_path`, `old_string`, `new_string`, `replace_all`   |
| `Read`  | `file_path`, `offset`, `limit`                           |

### PostToolUse Input

```json
{
  "session_id": "abc123",
  "transcript_path": "/path/to/session.jsonl",
  "cwd": "/path/to/project",
  "permission_mode": "default",
  "hook_event_name": "PostToolUse",
  "tool_name": "Write",
  "tool_input": {
    "file_path": "/path/to/file.txt",
    "content": "file content"
  },
  "tool_response": {
    "filePath": "/path/to/file.txt",
    "success": true
  },
  "tool_use_id": "toolu_01ABC123"
}
```

### Notification Input

```json
{
  "session_id": "abc123",
  "transcript_path": "/path/to/session.jsonl",
  "cwd": "/path/to/project",
  "permission_mode": "default",
  "hook_event_name": "Notification",
  "message": "Claude needs your permission to use Bash",
  "notification_type": "permission_prompt"
}
```

### UserPromptSubmit Input

```json
{
  "session_id": "abc123",
  "transcript_path": "/path/to/session.jsonl",
  "cwd": "/path/to/project",
  "permission_mode": "default",
  "hook_event_name": "UserPromptSubmit",
  "prompt": "Write a function to calculate factorial"
}
```

### Stop Input

```json
{
  "session_id": "abc123",
  "transcript_path": "/path/to/session.jsonl",
  "cwd": "/path/to/project",
  "permission_mode": "default",
  "hook_event_name": "Stop",
  "stop_hook_active": true
}
```

**`stop_hook_active`**: True when Claude is already continuing due to a stop hook. Check this to prevent infinite loops.

### SubagentStart Input

```json
{
  "session_id": "abc123",
  "transcript_path": "/path/to/session.jsonl",
  "cwd": "/path/to/project",
  "permission_mode": "default",
  "hook_event_name": "SubagentStart",
  "agent_id": "agent-abc123",
  "agent_type": "Explore"
}
```

**Fields**:

- `agent_id`: Unique identifier for the subagent
- `agent_type`: Agent name (built-in like "Bash", "Explore", "Plan", or custom agent names)

### SubagentStop Input

```json
{
  "session_id": "abc123",
  "transcript_path": "/path/to/session.jsonl",
  "cwd": "/path/to/project",
  "permission_mode": "default",
  "hook_event_name": "SubagentStop",
  "stop_hook_active": false,
  "agent_id": "def456",
  "agent_transcript_path": "/path/to/subagents/agent-def456.jsonl"
}
```

**Fields**:

- `agent_id`: Unique identifier for the subagent
- `agent_transcript_path`: Path to the subagent's own transcript in nested `subagents/` folder
- `stop_hook_active`: True when already continuing due to a stop hook

### PreCompact Input

```json
{
  "session_id": "abc123",
  "transcript_path": "/path/to/session.jsonl",
  "cwd": "/path/to/project",
  "permission_mode": "default",
  "hook_event_name": "PreCompact",
  "trigger": "manual",
  "custom_instructions": ""
}
```

### Setup Input

```json
{
  "session_id": "abc123",
  "transcript_path": "/path/to/session.jsonl",
  "cwd": "/path/to/project",
  "permission_mode": "default",
  "hook_event_name": "Setup",
  "trigger": "init"
}
```

**Fields**:

- `trigger`: Either `"init"` (from `--init` or `--init-only`) or `"maintenance"` (from `--maintenance`)
- Setup hooks have access to `CLAUDE_ENV_FILE` for persisting environment variables

### SessionStart Input

```json
{
  "session_id": "abc123",
  "transcript_path": "/path/to/session.jsonl",
  "cwd": "/path/to/project",
  "permission_mode": "default",
  "hook_event_name": "SessionStart",
  "source": "startup",
  "model": "claude-sonnet-4-20250514"
}
```

**Fields**:

- `source`: Indicates how session started: `"startup"` (new), `"resume"` (resumed), `"clear"` (after `/clear`), or `"compact"` (after compaction)
- `model`: Model identifier when available
- `agent_type`: (Optional) Present when Claude Code started with `claude --agent <name>`

### SessionEnd Input

```json
{
  "session_id": "abc123",
  "transcript_path": "/path/to/session.jsonl",
  "cwd": "/path/to/project",
  "permission_mode": "default",
  "hook_event_name": "SessionEnd",
  "reason": "exit"
}
```

**`reason` values**: `clear`, `logout`, `prompt_input_exit`, `other`

---

## Hook Output

### Exit Codes

| Code  | Behavior                                                         |
| ----- | ---------------------------------------------------------------- |
| 0     | Success. stdout processed (JSON or plain text)                   |
| 2     | Blocking error. stderr used as error message, fed back to Claude |
| Other | Non-blocking error. stderr shown in verbose mode (Ctrl+O)        |

**Important**: Claude Code does not see stdout if exit code is 0, except for `UserPromptSubmit` and `SessionStart` where stdout is added to context.

### Exit Code 2 Behavior Per Event

| Event                | Exit Code 2 Behavior                             |
| -------------------- | ------------------------------------------------ |
| `PreToolUse`         | Blocks tool call, shows stderr to Claude         |
| `PermissionRequest`  | Denies permission, shows stderr to Claude        |
| `PostToolUse`        | Shows stderr to Claude (tool already ran)        |
| `PostToolUseFailure` | Shows stderr to Claude (tool already failed)     |
| `Notification`       | Shows stderr to user only                        |
| `UserPromptSubmit`   | Blocks prompt, erases it, shows stderr to user   |
| `Stop`               | Blocks stoppage, shows stderr to Claude          |
| `SubagentStart`      | Shows stderr to user only                        |
| `SubagentStop`       | Blocks stoppage, shows stderr to Claude subagent |
| `PreCompact`         | Shows stderr to user only                        |
| `Setup`              | Shows stderr to user only                        |
| `SessionStart`       | Shows stderr to user only                        |
| `SessionEnd`         | Shows stderr to user only                        |

---

## JSON Output Control

**Important**: JSON output only processed with exit code 0. Exit code 2 uses stderr only.

### Common JSON Fields (All Events)

```json
{
  "continue": true,
  "stopReason": "Message shown when continue is false",
  "suppressOutput": false,
  "systemMessage": "Optional warning message shown to user"
}
```

| Field            | Type    | Effect                                           |
| ---------------- | ------- | ------------------------------------------------ |
| `continue`       | boolean | `false` stops Claude (takes precedence over all) |
| `stopReason`     | string  | Shown to user when `continue` is false           |
| `suppressOutput` | boolean | Hide stdout from transcript mode                 |
| `systemMessage`  | string  | Warning message shown to user                    |

**Precedence**: `continue: false` takes precedence over any `decision: "block"` output.

### PreToolUse JSON Output

```json
{
  "hookSpecificOutput": {
    "hookEventName": "PreToolUse",
    "permissionDecision": "allow",
    "permissionDecisionReason": "Auto-approved documentation file",
    "updatedInput": {
      "field_to_modify": "new value"
    },
    "additionalContext": "Current environment: production. Proceed with caution."
  }
}
```

| Field                      | Values                 | Effect                                     |
| -------------------------- | ---------------------- | ------------------------------------------ |
| `permissionDecision`       | `allow`, `deny`, `ask` | Controls tool execution                    |
| `permissionDecisionReason` | string                 | Shown to user (allow/ask) or Claude (deny) |
| `updatedInput`             | object                 | Modifies tool input before execution       |
| `additionalContext`        | string                 | Added to Claude's context                  |

**Note**: `decision` and `reason` fields are deprecated. Use `hookSpecificOutput.permissionDecision` and `hookSpecificOutput.permissionDecisionReason`. Deprecated `"approve"` and `"block"` map to `"allow"` and `"deny"`.

### PermissionRequest JSON Output

Allow with modified input:

```json
{
  "hookSpecificOutput": {
    "hookEventName": "PermissionRequest",
    "decision": {
      "behavior": "allow",
      "updatedInput": {
        "command": "npm run lint"
      }
    }
  }
}
```

Deny with message:

```json
{
  "hookSpecificOutput": {
    "hookEventName": "PermissionRequest",
    "decision": {
      "behavior": "deny",
      "message": "Command not allowed by policy",
      "interrupt": true
    }
  }
}
```

### PostToolUse JSON Output

```json
{
  "decision": "block",
  "reason": "Explanation for decision",
  "hookSpecificOutput": {
    "hookEventName": "PostToolUse",
    "additionalContext": "Additional information for Claude"
  }
}
```

- `"block"` automatically prompts Claude with `reason`
- `undefined` does nothing, `reason` is ignored

### UserPromptSubmit JSON Output

```json
{
  "decision": "block",
  "reason": "Prompt contains sensitive data",
  "hookSpecificOutput": {
    "hookEventName": "UserPromptSubmit",
    "additionalContext": "My additional context here"
  }
}
```

**Two ways to add context (exit code 0):**

1. **Plain text stdout**: Any non-JSON text is added as context
2. **JSON with `additionalContext`**: More structured control

**Blocking prompts:**

- `"decision": "block"` prevents prompt processing, erases prompt from context
- `"reason"` shown to user but not added to context

### Stop / SubagentStop JSON Output

```json
{
  "decision": "block",
  "reason": "Must be provided when Claude is blocked from stopping"
}
```

- `"block"` prevents Claude from stopping. Must populate `reason`
- `undefined` allows Claude to stop. `reason` is ignored

### Setup JSON Output

```json
{
  "hookSpecificOutput": {
    "hookEventName": "Setup",
    "additionalContext": "Repository initialized with custom configuration"
  }
}
```

**Note**: Multiple hooks' `additionalContext` values are concatenated. Setup hooks have access to `CLAUDE_ENV_FILE` for persisting environment variables.

### SessionStart JSON Output

```json
{
  "hookSpecificOutput": {
    "hookEventName": "SessionStart",
    "additionalContext": "Project context loaded successfully"
  }
}
```

**Note**: Multiple hooks' `additionalContext` values are concatenated.

---

## Sources

- [Hooks Reference](https://code.claude.com/docs/en/hooks.md) (accessed 2026-01-28)
- [Hooks Guide](https://code.claude.com/docs/en/hooks-guide.md)
- [Settings Reference](https://code.claude.com/docs/en/settings.md)
- [Plugin Components Reference](https://code.claude.com/docs/en/plugins-reference.md#hooks)
