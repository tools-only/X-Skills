# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/clay-pack/skills/clay-load-scale/SKILL.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/clay-pack/skills/clay-load-scale/SKILL.md` |
| **Category** | development |
| **Subcategory** | tools |
| **Tags** | development |
| **Created** | 2026-01-06 |
| **Updated** | 2026-01-06 |
| **File Hash** | `c778a28963124297...` |

## Description

export default function () {
  const response = http.post(
    'https://api.clay.com/v1/resource',
    JSON.stringify({ test: true }),
    {
      headers: {
        'ContentType': 'application/json',
        'Authorization': Bearer ${__ENV.CLAY_API_KEY},
      },
    }
  );

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/clay-pack/skills/clay-load-scale/SKILL.md)*
