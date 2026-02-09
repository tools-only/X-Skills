# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/openevidence-pack/skills/openevidence-common-errors/SKILL.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/openevidence-pack/skills/openevidence-common-errors/SKILL.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-01-26 |
| **Updated** | 2026-01-26 |
| **File Hash** | `bc3a0a08acc66e54...` |

## Description

static fromApiError(error: any): OpenEvidenceError {
    const statusCode = error.response?.status || 500;
    const message = error.response?.data?.message || error.message;
    const code = error.response?.data?.code || 'UNKNOWN_ERROR';

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/openevidence-pack/skills/openevidence-common-errors/SKILL.md)*
