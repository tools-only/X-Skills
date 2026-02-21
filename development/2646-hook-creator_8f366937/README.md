# Hook Creator

| Property | Value |
|----------|-------|
| **Name** | Hook Creator |
| **Repository** | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/plugin-creator/agents/hook-creator.md) (⭐ 18) |
| **Original Path** | `plugins/plugin-creator/agents/hook-creator.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-19 |
| **Updated** | 2026-02-19 |
| **File Hash** | `8f366937cfa346f5...` |

## Description

Creates Claude Code hook scripts for plugins — generates Node.js .cjs files, wires hooks.json, selects correct event and scope. Use when creating hooks, wiring PostToolUse or PreToolUse logic, enforcing validation on tool calls, or building SessionStart context injection. Trigger phrases — create a hook, add a hook to my plugin, build a PostToolUse hook, I need a hook that. <example>User asks to block rm -rf with a PreToolUse hook — hook-creator generates the .cjs script and wires hooks.json.</example> <example>User asks to inject project context on SessionStart — hook-creator builds the context injection script.</example> <example>User asks to run prettier after every Write — hook-creator wires a PostToolUse formatter hook.</example>

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/plugin-creator/agents/hook-creator.md)*
