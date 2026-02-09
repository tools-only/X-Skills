# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/lindy-pack/skills/lindy-observability/SKILL.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/lindy-pack/skills/lindy-observability/SKILL.md` |
| **Category** | development |
| **Subcategory** | tools |
| **Tags** | development |
| **Created** | 2026-01-09 |
| **Updated** | 2026-01-09 |
| **File Hash** | `ba404be0c15ca91c...` |

## Description

export const logger = pino({
  level: process.env.LOG_LEVEL || 'info',
  formatters: {
    level: (label) => ({ level: label }),
  },
  base: {
    service: 'lindyintegration',
    environment: process.env.NODE_ENV,
  },
});

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/lindy-pack/skills/lindy-observability/SKILL.md)*
