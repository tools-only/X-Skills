# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/lindy-pack/skills/lindy-multi-env-setup/SKILL.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/lindy-pack/skills/lindy-multi-env-setup/SKILL.md` |
| **Category** | development |
| **Subcategory** | devops |
| **Tags** | development |
| **Created** | 2026-01-09 |
| **Updated** | 2026-01-09 |
| **File Hash** | `3172ae0d50f55e84...` |

## Description

export function getLindyConfig(): LindyConfig {
  const env = process.env.NODE_ENV || 'development';
  return configs[env] || configs.development;
}

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/lindy-pack/skills/lindy-multi-env-setup/SKILL.md)*
