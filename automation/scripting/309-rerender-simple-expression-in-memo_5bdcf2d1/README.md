# Rerender Simple Expression In Memo

| Property | Value |
|----------|-------|
| **Name** | Rerender Simple Expression In Memo |
| **Repository** | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/vercel-react-best-practices/rules/rerender-simple-expression-in-memo.md) (⭐ 13) |
| **Original Path** | `skills/vercel-react-best-practices/rules/rerender-simple-expression-in-memo.md` |
| **Category** | automation |
| **Subcategory** | scripting |
| **Tags** | automation |
| **Created** | 2026-01-20 |
| **Updated** | 2026-01-20 |
| **File Hash** | `5bdcf2d1558ba520...` |

## Description

When an expression is simple (few logical or arithmetical operators) and has a primitive result type (boolean, number, string), do not wrap it in useMemo.
Calling useMemo and comparing hook dependencies may consume more resources than the expression itself.

**Tags:** `automation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/vercel-react-best-practices/rules/rerender-simple-expression-in-memo.md)*
