# Common Bug Patterns

| Property | Value |
|----------|-------|
| **Name** | Common Bug Patterns |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/debugging-wizard/references/common-patterns.md) (⭐ 216) |
| **Original Path** | `skills/debugging-wizard/references/common-patterns.md` |
| **Category** | commercial |
| **Subcategory** | ecommerce |
| **Tags** | commercial |
| **Created** | 2025-12-14 |
| **Updated** | 2026-01-29 |
| **File Hash** | `34e8957fc8dd4403...` |

## Description

typescript
// BUG: Race condition
let data;
fetchData().then(result => { data = result; });
console.log(data); // undefined!

**Tags:** `commercial`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/debugging-wizard/references/common-patterns.md)*
