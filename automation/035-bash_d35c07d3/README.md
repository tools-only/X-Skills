# Cloudflare D1 Kv

| Property | Value |
|----------|-------|
| **Name** | Cloudflare D1 Kv |
| **Repository** | [mrgoonie/claudekit-skills](https://raw.githubusercontent.com/mrgoonie/claudekit-skills/main/.claude/skills/devops/references/cloudflare-d1-kv.md) (⭐ 1.5k) |
| **Original Path** | `.claude/skills/devops/references/cloudflare-d1-kv.md` |
| **Category** | automation |
| **Subcategory** | scripting |
| **Tags** | automation |
| **Created** | 2025-11-06 |
| **Updated** | 2025-11-06 |
| **File Hash** | `d35c07d3bf5356a2...` |

## Description

typescript
// Query
const result = await env.DB.prepare(
  "SELECT  FROM users WHERE id = ?"
).bind(userId).first();

**Tags:** `automation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [mrgoonie/claudekit-skills](https://raw.githubusercontent.com/mrgoonie/claudekit-skills/main/.claude/skills/devops/references/cloudflare-d1-kv.md)*
