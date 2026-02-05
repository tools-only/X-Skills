# Rerender Use Ref Transient Values

| Property | Value |
|----------|-------|
| **Name** | Rerender Use Ref Transient Values |
| **Repository** | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/vercel-react-best-practices/rules/rerender-use-ref-transient-values.md) (⭐ 13) |
| **Original Path** | `skills/vercel-react-best-practices/rules/rerender-use-ref-transient-values.md` |
| **Category** | commercial |
| **Subcategory** | ecommerce |
| **Tags** | commercial |
| **Created** | 2026-01-29 |
| **Updated** | 2026-01-29 |
| **File Hash** | `f1a649af9d1d0b5c...` |

## Description

When a value changes frequently and you don't want a rerender on every update (e.g., mouse trackers, intervals, transient flags), store it in useRef instead of useState. Keep component state for UI; use refs for temporary DOMadjacent values. Updating a ref does not trigger a rerender.

**Tags:** `commercial`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/vercel-react-best-practices/rules/rerender-use-ref-transient-values.md)*
