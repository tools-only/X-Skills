# List Dir Tool Execution Error

| Property | Value |
|----------|-------|
| **Name** | List Dir Tool Execution Error |
| **Repository** | [alchemiststudiosDOTai/tunacode](https://raw.githubusercontent.com/alchemiststudiosDOTai/tunacode/master/.claude/debug_history/list-dir-tool-execution-error.md) (⭐ 112) |
| **Original Path** | `.claude/debug_history/list-dir-tool-execution-error.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-01-08 |
| **Updated** | 2026-01-08 |
| **File Hash** | `d1befeabae5d805d...` |

## Description

list_dir raised FileNotFoundError on nonexistent directories, which @base_tool wrapped as ToolExecutionError. Since ToolExecutionError is in NON_RETRYABLE_ERRORS, the agent halted instead of letting the LLM selfcorrect with a valid path.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [alchemiststudiosDOTai/tunacode](https://raw.githubusercontent.com/alchemiststudiosDOTai/tunacode/master/.claude/debug_history/list-dir-tool-execution-error.md)*
