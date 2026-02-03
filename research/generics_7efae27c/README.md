# Generics and Type Parameters

| Property | Value |
|----------|-------|
| **Name** | Generics and Type Parameters |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/golang-pro/references/generics.md) (⭐ 216) |
| **Original Path** | `skills/golang-pro/references/generics.md` |
| **Category** | research |
| **Subcategory** | academic |
| **Tags** | research |
| **Created** | 2025-12-15 |
| **Updated** | 2026-01-29 |
| **File Hash** | `7efae27c421afcd2...` |

## Description

// Generic function with type parameter
func Max[T constraints.Ordered](a, b T) T {
    if a > b {
        return a
    }
    return b
}

**Tags:** `research`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/golang-pro/references/generics.md)*
