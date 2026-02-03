# 2026 02 02 Tool Name Normalization

| Property | Value |
|----------|-------|
| **Name** | 2026 02 02 Tool Name Normalization |
| **Repository** | [alchemiststudiosDOTai/tunacode](https://raw.githubusercontent.com/alchemiststudiosDOTai/tunacode/master/.claude/delta/2026-02-02-tool-name-normalization.md) (⭐ 112) |
| **Original Path** | `.claude/delta/2026-02-02-tool-name-normalization.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-02 |
| **Updated** | 2026-02-02 |
| **File Hash** | `be6723f0224fa748...` |

## Description

Tool calls with leading/trailing whitespace in tool_name failed to dispatch (e.g., " glob"), resulting in unknowntool errors instead of retries. We now normalize tool names before registration and dispatch to avoid whitespaceinduced failures.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [alchemiststudiosDOTai/tunacode](https://raw.githubusercontent.com/alchemiststudiosDOTai/tunacode/master/.claude/delta/2026-02-02-tool-name-normalization.md)*
