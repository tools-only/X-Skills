# Api

| Property | Value |
|----------|-------|
| **Name** | Api |
| **Repository** | [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/do-storage/api.md) (⭐ 392) |
| **Original Path** | `packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/do-storage/api.md` |
| **Category** | communication |
| **Subcategory** | email |
| **Tags** | communication |
| **Created** | 2026-02-05 |
| **Updated** | 2026-02-05 |
| **File Hash** | `2038de52d380ed37...` |

## Description

typescript
const cursor = this.sql.exec('SELECT  FROM users WHERE email = ?', email)
for (let row of cursor) {
} // Objects: { id, name, email }
cursor.toArray()
cursor.one() // Single row (throws if != 1)
for (let row of cursor.raw()) {
} // Arrays: [1, "Alice", "..."]

**Tags:** `communication`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/do-storage/api.md)*
