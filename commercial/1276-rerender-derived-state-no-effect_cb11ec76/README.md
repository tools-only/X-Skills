# Rerender Derived State No Effect

| Property | Value |
|----------|-------|
| **Name** | Rerender Derived State No Effect |
| **Repository** | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/vercel-react-best-practices/rules/rerender-derived-state-no-effect.md) (⭐ 13) |
| **Original Path** | `skills/vercel-react-best-practices/rules/rerender-derived-state-no-effect.md` |
| **Category** | commercial |
| **Subcategory** | ecommerce |
| **Tags** | commercial |
| **Created** | 2026-01-29 |
| **Updated** | 2026-01-29 |
| **File Hash** | `cb11ec76f50aa7b6...` |

## Description

If a value can be computed from current props/state, do not store it in state or update it in an effect. Derive it during render to avoid extra renders and state drift. Do not set state in effects solely in response to prop changes; prefer derived values or keyed resets instead.

**Tags:** `commercial`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/vercel-react-best-practices/rules/rerender-derived-state-no-effect.md)*
