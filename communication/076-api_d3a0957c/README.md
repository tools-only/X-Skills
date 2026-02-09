# Api

| Property | Value |
|----------|-------|
| **Name** | Api |
| **Repository** | [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/do-storage/api.md) (🔥 19.8k) |
| **Original Path** | `cli-tool/components/skills/development/cloudflare-deploy/references/do-storage/api.md` |
| **Category** | communication |
| **Subcategory** | email |
| **Tags** | communication |
| **Created** | 2026-02-08 |
| **Updated** | 2026-02-08 |
| **File Hash** | `d3a0957cca017ae8...` |

## Description

typescript
const cursor = this.sql.exec('SELECT  FROM users WHERE email = ?', email);
for (let row of cursor) {} // Objects: { id, name, email }
cursor.toArray(); cursor.one(); // Single row (throws if != 1)
for (let row of cursor.raw()) {} // Arrays: [1, "Alice", "..."]

**Tags:** `communication`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/do-storage/api.md)*
