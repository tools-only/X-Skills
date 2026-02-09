# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/gamma-pack/skills/gamma-data-handling/SKILL.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/gamma-pack/skills/gamma-data-handling/SKILL.md` |
| **Category** | content-creation |
| **Subcategory** | writing |
| **Tags** | content creation |
| **Created** | 2026-01-09 |
| **Updated** | 2026-01-09 |
| **File Hash** | `a912653e85e32ad6...` |

## Description

async function checkConsent(userId: string, purpose: string): Promise<boolean> {
  const consent = await db.consents.findUnique({
    where: { userId },
  });

**Tags:** `content creation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/gamma-pack/skills/gamma-data-handling/SKILL.md)*
