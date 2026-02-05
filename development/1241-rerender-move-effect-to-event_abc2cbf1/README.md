# Rerender Move Effect To Event

| Property | Value |
|----------|-------|
| **Name** | Rerender Move Effect To Event |
| **Repository** | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/vercel-react-best-practices/rules/rerender-move-effect-to-event.md) (⭐ 13) |
| **Original Path** | `skills/vercel-react-best-practices/rules/rerender-move-effect-to-event.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-01-29 |
| **Updated** | 2026-01-29 |
| **File Hash** | `abc2cbf167bee056...` |

## Description

If a side effect is triggered by a specific user action (submit, click, drag), run it in that event handler. Do not model the action as state + effect; it makes effects rerun on unrelated changes and can duplicate the action.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/vercel-react-best-practices/rules/rerender-move-effect-to-event.md)*
