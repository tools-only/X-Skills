# Hook Patterns for Kaizen Improvements

Reference for generating Claude Code hooks from analysis findings. Each pattern maps an anti-pattern to a hook configuration that prevents or corrects it.

## PreToolUse Hooks — Deny and Redirect

### Tool Misuse Prevention

When analysis identifies Bash used for file operations, generate a PreToolUse hook that denies the Bash call and tells Claude to use the correct tool.

```json
{
  "PreToolUse": [{
    "matcher": "Bash",
    "hooks": [{
      "type": "command",
      "command": "node ${CLAUDE_PLUGIN_ROOT}/hooks/scripts/check-tool-misuse.js"
    }]
  }]
}
```

Hook script pattern — check `tool_input.command` for violations:

```javascript
#!/usr/bin/env node
const input = JSON.parse(require('fs').readFileSync('/dev/stdin', 'utf8'));
const cmd = input.tool_input?.command || '';

const violations = [
  { pattern: /^\s*ls\b/, redirect: 'Glob', msg: 'Use Glob tool instead of ls' },
  { pattern: /\bgrep\b(?!.*\|.*grep)/, redirect: 'Grep', msg: 'Use Grep tool instead of grep' },
  { pattern: /\bcat\b\s+\S+\.\w+/, redirect: 'Read', msg: 'Use Read tool instead of cat' },
  { pattern: /\bfind\b\s+\S+\s+-name\b/, redirect: 'Glob', msg: 'Use Glob tool instead of find' },
  { pattern: /\bhead\b\s+-\d+/, redirect: 'Read', msg: 'Use Read tool with limit instead of head' },
  { pattern: /\btail\b\s+-\d+/, redirect: 'Read', msg: 'Use Read tool with offset instead of tail' },
];

for (const v of violations) {
  if (v.pattern.test(cmd)) {
    console.log(JSON.stringify({
      hookSpecificOutput: {
        hookEventName: 'PreToolUse',
        permissionDecision: 'deny',
        permissionDecisionReason: v.msg,
        additionalContext: `Redirect to ${v.redirect} tool. Original command: ${cmd}`
      }
    }));
    process.exit(0);
  }
}
process.exit(0);
```

### Edit-Before-Read Prevention

Deny Edit calls when the target file has not been Read in the current session.

**Caveat:** A `prompt` type hook relies on the model's session memory, which may be unreliable after context compaction. For deterministic enforcement, use a `command` type hook with state tracking (e.g., a script that records Read calls and checks against them).

```json
{
  "PreToolUse": [{
    "matcher": "Edit",
    "hooks": [{
      "type": "prompt",
      "prompt": "Check if the file in tool_input.file_path has been Read in this session. If not, deny with reason 'Read the file before editing it.'"
    }]
  }]
}
```

## SubagentStart Hooks — Inject Context

### Research-to-Files Enforcement

When analysis shows subagents sending research as messages instead of writing to files:

```json
{
  "SubagentStart": [{
    "hooks": [{
      "type": "prompt",
      "prompt": "Inject this context into the subagent: 'CRITICAL: Write all research findings to files, not messages. Messages are lost during context compaction. Write to .claude/ or .planning/ directory.'"
    }]
  }]
}
```

### Delegation Pattern Enforcement

When analysis shows orchestrators dictating solutions instead of delegating outcomes:

```json
{
  "SubagentStart": [{
    "hooks": [{
      "type": "prompt",
      "prompt": "Check the subagent's task prompt. If it contains specific implementation instructions (exact code, specific variable names, line-by-line changes), deny with: 'Delegation prompts must describe outcomes, not solutions. Describe the problem and desired outcome, let the agent decide the approach.'"
    }]
  }]
}
```

## SubagentStop Hooks — Quality Gates

### Output Validation

When analysis shows subagents producing incomplete or low-quality output:

```json
{
  "SubagentStop": [{
    "hooks": [{
      "type": "prompt",
      "prompt": "Validate the subagent's output. Check: (1) Did it write findings to files? (2) Did it cite sources? (3) Does the output match the task description? Return ok:false with specific gaps if validation fails."
    }]
  }]
}
```

## Stop Hooks — Completion Verification

### Commit-Before-Stop

When analysis shows work being lost due to missing commits:

```json
{
  "Stop": [{
    "hooks": [{
      "type": "prompt",
      "prompt": "Check if there are uncommitted changes related to the completed task. If so, suggest: 'Consider committing your changes before stopping. Uncommitted work can be lost during context compaction or session restart.'"
    }]
  }]
}
```

## Hook Generation Guidelines

When translating findings into hooks:

1. **One anti-pattern per hook** — keep hooks focused and testable
2. **Prefer `command` type for deterministic checks** — regex matching, file existence
3. **Use `prompt` type for semantic checks** — understanding intent, quality assessment
4. **Include `matcher` to scope narrowly** — only trigger on relevant tools
5. **Test with `claude --debug`** — verify hook fires and produces expected output
6. **Draft first, install later** — write to `.planning/kaizen/hooks/` for review before `--install`

SOURCE: Claude Code hooks documentation (code.claude.com/docs/en/hooks.md, code.claude.com/docs/en/hooks-guide.md), accessed 2026-02-17
