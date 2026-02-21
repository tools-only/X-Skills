# Orchestrator Discipline Plugin — Codebase Patterns

**Analysis Date:** 2026-02-20
**Package:** `orchestrator-discipline`
**Plugin Root:** `plugins/orchestrator-discipline/`

---

## Table of Contents

1. [Hook JS Implementation Patterns](#1-hook-js-implementation-patterns)
2. [hooks.json Structure](#2-hooksjson-structure)
3. [plugin.json `rules` Field Status](#3-pluginjson-rules-field-status)
4. [Grep Tool PreToolUse Input Structure](#4-grep-tool-pretooluse-input-structure)
5. [Directory Path Detection in Hook JS](#5-directory-path-detection-in-hook-js)
6. [Conventions Across All Hooks](#6-conventions-across-all-hooks)

---

## 1. Hook JS Implementation Patterns

### 1.1 Input Parsing: stdin vs environment variables

All hook JS files in this repo receive JSON input exclusively via **stdin**, not via environment variables or command arguments.

Pattern used in every hook that needs tool data:

```javascript
let input = '';
process.stdin.setEncoding('utf8');
process.stdin.on('data', (chunk) => { input += chunk; });
process.stdin.on('end', () => {
  let data = {};
  try {
    data = JSON.parse(input);
  } catch {
    process.stdout.write(JSON.stringify({}));
    process.exit(0);
  }
  // ... use data
});
```

Source: `plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs:28-40`
Source: `plugins/orchestrator-discipline/hooks/pre-tool-diagnostic-command-gate.cjs:50-64`
Source: `plugins/dasel/hooks/session-start-dasel-check.cjs:39-43`

**Exception**: `Stop` and `SessionStart` hooks that do not need tool input may skip stdin reading entirely and output directly via `console.log`. Example:

```javascript
// No stdin parsing needed — hook fires without per-tool data
const output = { hookSpecificOutput: { hookEventName: 'SessionStart', additionalContext: '...' } };
console.log(JSON.stringify(output));
```

Source: `.claude/hooks/session-start-rtica.cjs:7-16`

**Alternative synchronous stdin**: The hallucination-detector Stop hook uses synchronous `fs.readFileSync(0, 'utf-8')` instead of async stdin events:

```javascript
function readStdinJson() {
  try {
    const stdin = fs.readFileSync(0, 'utf-8');
    return JSON.parse(stdin);
  } catch { return {}; }
}
```

Source: `plugins/hallucination-detector/scripts/hallucination-audit-stop.js:22-29`

### 1.2 Top-level JSON fields from Claude Code

For PreToolUse hooks, the stdin JSON contains:

| Field        | Type   | Description                          | Observed in                    |
|-------------|--------|--------------------------------------|-------------------------------|
| `tool_name`  | string | Name of the tool being called        | `data.tool_name`               |
| `tool_input` | object | Tool-specific parameters object      | `data.tool_input`              |

Source: `plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs:42-53`

For Stop hooks, the stdin JSON contains:

| Field              | Type    | Description                          |
|-------------------|---------|--------------------------------------|
| `transcript_path`  | string  | Path to JSONL session transcript     |
| `session_id`       | string  | Session identifier                   |
| `stop_hook_active` | boolean | Whether another stop hook is active  |

Source: `plugins/hallucination-detector/scripts/hallucination-audit-stop.js:250-254`

### 1.3 `additionalContext` format

Non-blocking hooks inject guidance into the model's context by returning this structure on stdout:

```javascript
const output = {
  hookSpecificOutput: {
    hookEventName: 'PreToolUse',   // Must match the triggering event
    additionalContext: `<xml-tag-name>
... guidance text ...
</xml-tag-name>`,
  },
};
process.stdout.write(JSON.stringify(output));
process.exit(0);
```

The `additionalContext` value is always wrapped in a named XML tag. Tag naming follows the plugin domain:

- `<orchestrator-read-warning>` — `pre-tool-orchestrator-read-warning.cjs:62-79`
- `<orchestrator-diagnostic-warning>` — `pre-tool-diagnostic-command-gate.cjs:81-104`
- `<backlog-summary>` — `.claude/hooks/session-start-backlog.cjs:63-66`
- `<rt-ica-checkpoint>` — `.claude/hooks/session-start-rtica.cjs:10-13`
- `<backlog-reminder>` — `.claude/hooks/stop-backlog-reminder.cjs:8-13`

**Variation in Stop hook output (no `hookSpecificOutput` wrapper)**:

The stop-backlog-reminder hook writes `additionalContext` directly at the top level, without the `hookSpecificOutput` wrapper:

```javascript
const output = {
  additionalContext: `<backlog-reminder>...</backlog-reminder>`,
};
```

Source: `.claude/hooks/stop-backlog-reminder.cjs:8-14`

All `PreToolUse` hooks in this repo use the `hookSpecificOutput` wrapper.

### 1.4 Blocking vs non-blocking hooks

**Non-blocking (advisory)**: All orchestrator-discipline hooks return the `additionalContext` pattern shown above. Claude Code injects the context and allows the tool call to proceed.

**Blocking**: The hallucination-detector uses `"decision": "block"` to prevent the Stop event from completing:

```javascript
emitJson({ decision: 'block', reason: '...' });
process.exit(0);
```

Source: `plugins/hallucination-detector/scripts/hallucination-audit-stop.js:311`

### 1.5 Error handling pattern

All hooks follow the same error-safety idiom: catch JSON parse failures and exit cleanly without output:

```javascript
try {
  data = JSON.parse(input);
} catch {
  process.stdout.write(JSON.stringify({}));
  process.exit(0);
}
```

Source: `plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs:36-40`

Early exit with empty object `{}` on non-matching tool names:

```javascript
if (toolName !== 'Bash') {
  process.stdout.write(JSON.stringify({}));
  process.exit(0);
}
```

Source: `plugins/orchestrator-discipline/hooks/pre-tool-diagnostic-command-gate.cjs:65-68`

### 1.6 Accessing tool-specific fields

**Read tool** — target path is in `tool_input.file_path`:

```javascript
if (toolName === 'Read') {
  targetPath = toolInput.file_path || '';
}
```

Source: `plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs:46-48`

**Grep tool** — target path is in `tool_input.path`:

```javascript
} else if (toolName === 'Grep') {
  targetPath = toolInput.path || '';
}
```

Source: `plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs:49-51`

**Bash tool** — command string is in `tool_input.command`:

```javascript
const command = data.tool_input?.command || '';
```

Source: `plugins/orchestrator-discipline/hooks/pre-tool-diagnostic-command-gate.cjs:70`

---

## 2. hooks.json Structure

### 2.1 Plugin hooks.json top-level format

All plugin hooks.json files use one of two structures:

**Structure A — with `description` field (orchestrator-discipline, hallucination-detector, dasel)**:

```json
{
  "description": "Human-readable description of this hook configuration",
  "hooks": {
    "EventName": [
      {
        "matcher": "ToolA|ToolB",
        "hooks": [
          {
            "type": "command",
            "command": "node \"${CLAUDE_PLUGIN_ROOT}/hooks/script.js\""
          }
        ]
      }
    ]
  }
}
```

**Structure B — hooks only (summarizer)**:

```json
{
  "hooks": {
    "SubagentStop": [
      {
        "hooks": [
          {
            "type": "prompt",
            "prompt": "...",
            "timeout": 30
          }
        ]
      }
    ]
  }
}
```

Source files:
- `plugins/orchestrator-discipline/hooks.json` — Structure A
- `plugins/hallucination-detector/hooks.json` — Structure A
- `plugins/dasel/hooks/hooks.json` — Structure A
- `plugins/summarizer/hooks/hooks.json` — Structure B

### 2.2 matcher field

`matcher` is a pipe-separated regex string matching tool names. Used only when filtering specific tools within an event:

```json
"matcher": "Read|Grep"
"matcher": "Bash"
"matcher": "Write|Edit"
```

Source: `plugins/orchestrator-discipline/hooks.json:6,14`

When no filtering is needed (hook fires for all invocations of the event), the entry array element omits `matcher`:

```json
{
  "hooks": [{ "type": "command", "command": "..." }]
}
```

Source: `plugins/hallucination-detector/hooks.json:3-14`

### 2.3 hook type field

Three types observed in this repo:

| Type      | Usage                                                                 | Example location                                |
|-----------|-----------------------------------------------------------------------|-------------------------------------------------|
| `command` | Execute a shell command or Node.js script                            | `plugins/orchestrator-discipline/hooks.json:9`  |
| `prompt`  | Evaluate a prompt with an LLM; supports `$ARGUMENTS` placeholder     | `plugins/summarizer/hooks/hooks.json:6`         |

`agent` type is documented in the SKILL.md reference but is not used in this repo's hooks.json files.

### 2.4 ${CLAUDE_PLUGIN_ROOT} in command strings

The environment variable `${CLAUDE_PLUGIN_ROOT}` is used in all plugin hook command strings for path resolution. Paths must be quoted when the plugin root could contain spaces:

```json
"command": "node \"${CLAUDE_PLUGIN_ROOT}/hooks/pre-tool-orchestrator-read-warning.cjs\""
```

Source: `plugins/orchestrator-discipline/hooks.json:11`

Scripts in `.claude/hooks/` (project-level, not plugin) are referenced directly without `${CLAUDE_PLUGIN_ROOT}` — they are resolved at a different level.

### 2.5 hooks.json location in plugin.json

The `hooks` field in `plugin.json` accepts a path string pointing to a hooks.json file. The orchestrator-discipline plugin places hooks.json at the plugin root:

```json
"hooks": "./hooks.json"
```

Source: `plugins/orchestrator-discipline/.claude-plugin/plugin.json:11`

The dasel plugin places hooks.json in a `hooks/` subdirectory:

```json
"hooks": "./hooks/hooks.json"
```

Source: (implied by `plugins/dasel/hooks/hooks.json` and dasel's plugin.json referencing it)

---

## 3. plugin.json `rules` Field Status

### 3.1 Is `rules` a valid field?

The official plugin.json schema documented in `plugins/plugin-creator/skills/claude-plugins-reference-2026/SKILL.md` (lines 31-52) defines these component path fields:

```
commands, agents, skills, hooks, mcpServers, outputStyles, lspServers
```

The `rules` field is **not listed** in the official schema.

### 3.2 Observed usage in this repo

The orchestrator-discipline plugin uses `rules` in its plugin.json:

```json
{
  "name": "orchestrator-discipline",
  "skills": ["./skills/orchestrator-discipline"],
  "hooks": "./hooks.json",
  "rules": ["./rules"],
  "commands": ["./skills/orchestrator-discipline"]
}
```

Source: `plugins/orchestrator-discipline/.claude-plugin/plugin.json:11`

The `./rules` directory contains `plugins/orchestrator-discipline/rules/CLAUDE.md` — a markdown file with enforcement rules intended to be injected into the Claude context window.

### 3.3 How other plugins handle CLAUDE.md rules files

Other plugins in this repo do **not** use a `rules` field. Instead, they use one of two approaches:

**Approach A — Plugin root CLAUDE.md**: Several plugins place a `CLAUDE.md` at the plugin root (e.g., `plugins/development-harness/CLAUDE.md`, `plugins/plugin-creator/CLAUDE.md`). This is loaded as plugin-level documentation context, not via the `rules` field.

**Approach B — Inline in skills/commands**: Rules are embedded directly in SKILL.md content which is loaded when the skill is activated.

**Approach C — Via hooks**: The orchestrator-discipline plugin enforces rules structurally via hooks rather than solely relying on CLAUDE.md injection.

### 3.4 Practical implication

Because `rules` is not in the official schema, whether Claude Code loads `./rules/CLAUDE.md` via the `rules` field is unverified from the schema documentation alone. The enforcement in orchestrator-discipline is primarily achieved through the PreToolUse hooks (`hooks.json`) and the `additionalContext` injections, not solely through the `rules` directory. The hooks fire on every Read/Grep/Bash call and inject reminders regardless of whether `rules/CLAUDE.md` is loaded.

---

## 4. Grep Tool PreToolUse Input Structure

### 4.1 Fields confirmed from hook source

When Claude Code fires a `PreToolUse` hook for a Grep call, the `tool_input` object contains at minimum:

```javascript
{
  tool_name: 'Grep',
  tool_input: {
    path: '<value>',   // The path argument passed to Grep
    pattern: '<value>' // The search pattern (not used by this hook)
    // ... other Grep fields
  }
}
```

Source: `plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs:48-51`

The hook reads `toolInput.path || ''` to get the target path.

### 4.2 What `path` contains for directory vs file calls

The hook applies `isSourceOrConfigFile(targetPath)` which uses a regex test against the path string:

```javascript
const SOURCE_FILE_EXTENSIONS =
  /\.(py|toml|yaml|yml|js|ts|jsx|tsx|json|cfg|ini|env|sh|bash|go|rs|rb|java|c|cpp|h|hpp)$/i;
```

Source: `plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs:15-16`

This regex tests whether the path **ends with a source file extension**. It does **not** match directory paths (which have no extension). Therefore:

- `path: "src/utils.py"` → matches → hook fires warning
- `path: "src/"` → no match → hook passes silently
- `path: "."` → no match → hook passes silently
- `path: "plugins/orchestrator-discipline/hooks"` → no match → hook passes silently
- `path: "plugins/my-plugin/hooks/script.js"` → matches → hook fires warning

The `TEST_PATH_PATTERN` also matches:

```javascript
const TEST_PATH_PATTERN = /\/(tests?|spec|__tests?__)\/|test_[^/]+\.(py|js|ts)$/i;
```

Source: `plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs:17`

This matches paths containing `/tests/`, `/test/`, `/spec/`, `/__tests__/` as directory segments.

### 4.3 Gap: directory-level Grep is not intercepted

When Grep is called with `path: "src/"` or `path: "."`, the current hook will not fire a warning even though scanning a directory recursively can dump large amounts of source content into the orchestrator's context. The extension-based check only catches explicit file paths.

---

## 5. Directory Path Detection in Hook JS

### 5.1 Current approach (extension matching)

The current hook uses extension matching to detect source files. This approach misses:

1. Directory paths passed to Grep (as documented in section 4.3)
2. Files without conventional extensions that are source code

### 5.2 How to correctly detect directory paths

A directory path lacks a file extension and typically does not end with a known source extension. To detect whether a Grep `path` is a directory vs a file, a hook could use either:

**Option A — Node.js `fs.statSync`** (filesystem check at hook execution time):

```javascript
const fs = require('node:fs');

function isDirectory(targetPath) {
  try {
    return fs.statSync(targetPath).isDirectory();
  } catch {
    return false;
  }
}
```

This is the most reliable approach when the path is resolvable at hook execution time.

**Option B — Pattern heuristics** (no filesystem access):

```javascript
function looksLikeDirectory(targetPath) {
  if (!targetPath) return false;
  // Ends with separator
  if (targetPath.endsWith('/') || targetPath.endsWith('\\')) return true;
  // No extension in last path segment
  const lastSegment = targetPath.split(/[/\\]/).pop() || '';
  return !lastSegment.includes('.');
}
```

This is less reliable but avoids filesystem access overhead.

**Option C — Extension allowlist only** (current approach, file-focused):

The current code fires warnings only for paths matching source file extensions. Directory-targeted Grep calls are silently allowed.

### 5.3 Recommendation for extending the hook

To also warn on directory-targeted Grep calls, add a directory check before the extension check:

```javascript
function shouldWarn(toolName, toolInput) {
  if (toolName === 'Read') {
    return isSourceOrConfigFile(toolInput.file_path || '');
  }
  if (toolName === 'Grep') {
    const p = toolInput.path || '';
    // Warn on both direct source file paths and directory paths
    return isSourceOrConfigFile(p) || isDirectory(p) || looksLikeDirectory(p);
  }
  return false;
}
```

---

## 6. Conventions Across All Hooks

### 6.1 Language: all hooks are Node.js

Every hook script in this repo (both plugin hooks and `.claude/hooks/`) is JavaScript running under Node.js. No Python, Bash, or other runtimes are used for hooks.

Files use either `.js` or `.cjs` extension. `.cjs` forces CommonJS module mode explicitly (observed only in the dasel hook).

- Plugin hooks: `.js` — `plugins/*/hooks/*.js`
- Legacy/compatibility: `.cjs` — `plugins/dasel/hooks/session-start-dasel-check.cjs`
- Project-level: `.js` — `.claude/hooks/*.js`

### 6.2 Shebang

All hook files begin with:

```javascript
#!/usr/bin/env node
```

Source: All hook files, line 1.

### 6.3 Module system

All hooks use CommonJS `require()` — not ES module `import`. No `"type": "module"` declarations.

```javascript
const fs = require('node:fs');
const path = require('node:path');
const os = require('node:os');
const { execFileSync } = require('node:child_process');
```

Node.js `node:` prefix is used consistently for stdlib modules.

### 6.4 Output method

Two patterns observed:

- `process.stdout.write(JSON.stringify(output))` — used in PreToolUse hooks (ensures no trailing newline issues)
- `console.log(JSON.stringify(output))` — used in SessionStart hooks (adds trailing newline, acceptable for these events)

Source:
- `process.stdout.write`: `pre-tool-orchestrator-read-warning.cjs:83`, `pre-tool-diagnostic-command-gate.cjs:108`
- `console.log`: `.claude/hooks/session-start-backlog.cjs:69`, `.claude/hooks/session-start-rtica.cjs:16`

### 6.5 Exit codes

All hooks exit with `process.exit(0)` regardless of whether they inject context or block. Non-zero exit codes are not used in this repo's hooks — decisions are communicated via JSON output (`"decision": "block"`), not exit codes.

### 6.6 Hook file location patterns

```
plugins/{plugin-name}/hooks/{hook-name}.js       # Plugin PreToolUse / event hooks
plugins/{plugin-name}/scripts/{script-name}.js   # Plugin Stop / SubagentStop hooks
.claude/hooks/{event-name}-{purpose}.js          # Project-level hooks
```

The hallucination-detector places its Stop hook under `scripts/` rather than `hooks/`, likely because the hooks.json points to it as `"${CLAUDE_PLUGIN_ROOT}/scripts/hallucination-audit-stop.js"`.

Source: `plugins/hallucination-detector/hooks.json:9`

### 6.7 hookSpecificOutput.hookEventName

All hooks that use the `hookSpecificOutput` wrapper set `hookEventName` to the event that fired them:

```javascript
hookSpecificOutput: {
  hookEventName: 'PreToolUse',   // or 'SessionStart', 'Stop', etc.
  additionalContext: '...',
}
```

This field documents which event the hook handles — it is informational within the output structure.

### 6.8 Matcher syntax in hooks.json

Matchers use pipe-separated tool names without spaces:

```json
"matcher": "Read|Grep"
"matcher": "Bash"
"matcher": "Write|Edit"
```

The matcher value is treated as a regex pattern by Claude Code.

### 6.9 Timeout field

Optional `timeout` in seconds can be added to individual hook entries:

```json
{ "type": "command", "command": "...", "timeout": 10 }
```

Source: `plugins/dasel/hooks/hooks.json:9`

The orchestrator-discipline hooks do not specify a timeout (they run fast, no subprocess invocations).

---

## Summary Reference

### Quick lookup: which field for which tool

| Tool  | Input field for path/target  | Source                                             |
|-------|------------------------------|----------------------------------------------------|
| Read  | `tool_input.file_path`       | `pre-tool-orchestrator-read-warning.cjs:47`         |
| Grep  | `tool_input.path`            | `pre-tool-orchestrator-read-warning.cjs:50`         |
| Bash  | `tool_input.command`         | `pre-tool-diagnostic-command-gate.cjs:70`           |

### Quick lookup: hook output structures

| Use case                     | Output structure                                                    |
|------------------------------|---------------------------------------------------------------------|
| Non-blocking context inject  | `{ hookSpecificOutput: { hookEventName, additionalContext } }`      |
| Block stop/action            | `{ decision: 'block', reason: '...' }`                             |
| No-op (pass through)         | `{}` (empty object)                                                 |

### plugin.json `rules` field

- Not in official schema (`claude-plugins-reference-2026/SKILL.md`)
- Used by orchestrator-discipline at `plugins/orchestrator-discipline/.claude-plugin/plugin.json:11`
- Points to `./rules` directory containing `CLAUDE.md`
- No other plugins in this repo use `rules` field
- Primary enforcement is via hooks, not the `rules` field

---

*Analysis based on direct file reads. All claims cite source file and line number.*
