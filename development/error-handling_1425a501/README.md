# Alice

| Property | Value |
|----------|-------|
| **Name** | Alice |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/rust-engineer/references/error-handling.md) (⭐ 216) |
| **Original Path** | `skills/rust-engineer/references/error-handling.md` |
| **Category** | data-analysis |
| **Subcategory** | processing |
| **Tags** | data analysis |
| **Created** | 2025-12-15 |
| **Updated** | 2026-01-29 |
| **File Hash** | `1425a501e3201760...` |

## Description

rust
// Result: operation that can fail
fn divide(a: f64, b: f64) > Result<f64, String> {
    if b == 0.0 {
        Err("Division by zero".to_string())
    } else {
        Ok(a / b)
    }
}

**Tags:** `data analysis`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/rust-engineer/references/error-handling.md)*
