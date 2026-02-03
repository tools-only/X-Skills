# XSS & CSRF Prevention

| Property | Value |
|----------|-------|
| **Name** | XSS & CSRF Prevention |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/secure-code-guardian/references/xss-csrf.md) (⭐ 216) |
| **Original Path** | `skills/secure-code-guardian/references/xss-csrf.md` |
| **Category** | data-analysis |
| **Subcategory** | processing |
| **Tags** | data analysis |
| **Created** | 2025-12-14 |
| **Updated** | 2026-01-29 |
| **File Hash** | `4d9976faed6ff05c...` |

## Description

typescript
// React automatically escapes by default
function SafeComponent({ userInput }: { userInput: string }) {
  return <div>{userInput}</div>; // Safe  autoescaped
}

**Tags:** `data analysis`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/secure-code-guardian/references/xss-csrf.md)*
