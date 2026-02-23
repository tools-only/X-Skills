# 2026 02 21 Minimax Tool Call Bugs

| Property | Value |
|----------|-------|
| **Name** | 2026 02 21 Minimax Tool Call Bugs |
| **Repository** | [alchemiststudiosDOTai/tunacode](https://raw.githubusercontent.com/alchemiststudiosDOTai/tunacode/master/.claude/delta/2026-02-21_minimax-tool-call-bugs.md) (⭐ 113) |
| **Original Path** | `.claude/delta/2026-02-21_minimax-tool-call-bugs.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | minimax, tool-call, bug |
| **Created** | 2026-02-21 |
| **Updated** | 2026-02-21 |
| **File Hash** | `3aa542efaf1a3482...` |

## Description

MiniMax M2.5 returns tool calls with empty id and sometimes empty name fields. Our Python code rejects empty strings in three places. The tinyagent Rust layer (alchemyllm) handles IDs correctly via ToolCallId with Into trait  this is NOT a tinyagent bug.

**Tags:** `minimax` `tool-call` `bug`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [alchemiststudiosDOTai/tunacode](https://raw.githubusercontent.com/alchemiststudiosDOTai/tunacode/master/.claude/delta/2026-02-21_minimax-tool-call-bugs.md)*
