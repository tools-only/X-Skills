# 2026 01 17 Dangling Tool Calls

| Property | Value |
|----------|-------|
| **Name** | 2026 01 17 Dangling Tool Calls |
| **Repository** | [alchemiststudiosDOTai/tunacode](https://raw.githubusercontent.com/alchemiststudiosDOTai/tunacode/master/.claude/delta/2026-01-17-dangling-tool-calls.md) (⭐ 112) |
| **Original Path** | `.claude/delta/2026-01-17-dangling-tool-calls.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-01-17 |
| **Updated** | 2026-01-17 |
| **File Hash** | `e44184fc8ea9ce3c...` |

## Description

When a user aborted midtoolcall (Ctrl+C or tool denial), the conversation history was left with a ModelResponse containing tool calls but no corresponding ToolReturn messages. The next API request failed because the provider expected tool returns for pending calls.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [alchemiststudiosDOTai/tunacode](https://raw.githubusercontent.com/alchemiststudiosDOTai/tunacode/master/.claude/delta/2026-01-17-dangling-tool-calls.md)*
