# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/clerk-pack/skills/clerk-migration-deep-dive/SKILL.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/clerk-pack/skills/clerk-migration-deep-dive/SKILL.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-01-09 |
| **Updated** | 2026-01-09 |
| **File Hash** | `6dff778f4a2ce55c...` |

## Description

interface ClerkImportUser {
  external_id: string
  email_addresses: Array<{
    email_address: string
    verified: boolean
  }>
  first_name?: string
  last_name?: string
  image_url?: string
  created_at: string
  public_metadata?: Record<string, any>
}

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/clerk-pack/skills/clerk-migration-deep-dive/SKILL.md)*
