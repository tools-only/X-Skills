# State Decouple Implementation

| Property | Value |
|----------|-------|
| **Name** | State Decouple Implementation |
| **Repository** | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/vercel-composition-patterns/rules/state-decouple-implementation.md) (⭐ 13) |
| **Original Path** | `skills/vercel-composition-patterns/rules/state-decouple-implementation.md` |
| **Category** | automation |
| **Subcategory** | scripting |
| **Tags** | automation |
| **Created** | 2026-01-29 |
| **Updated** | 2026-01-29 |
| **File Hash** | `76a8c1eb29a70d89...` |

## Description

The provider component should be the only place that knows how state is managed.
UI components consume the context interface—they don't know if state comes from
useState, Zustand, or a server sync.

**Tags:** `automation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/vercel-composition-patterns/rules/state-decouple-implementation.md)*
