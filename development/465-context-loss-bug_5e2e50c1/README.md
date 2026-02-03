# Context Loss Bug

| Property | Value |
|----------|-------|
| **Name** | Context Loss Bug |
| **Repository** | [alchemiststudiosDOTai/tunacode](https://raw.githubusercontent.com/alchemiststudiosDOTai/tunacode/master/.claude/debug_history/context-loss-bug.md) (⭐ 112) |
| **Original Path** | `.claude/debug_history/context-loss-bug.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-01-23 |
| **Updated** | 2026-01-23 |
| **File Hash** | `5e2e50c1b8b31f77...` |

## Description

Messages are only saved after successful completion of the agent loop. Any error, abort, or cancellation causes ALL messages to be lost, even if the conversation ran for many iterations.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [alchemiststudiosDOTai/tunacode](https://raw.githubusercontent.com/alchemiststudiosDOTai/tunacode/master/.claude/debug_history/context-loss-bug.md)*
