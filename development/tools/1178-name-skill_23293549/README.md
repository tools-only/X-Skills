# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/windsurf-pack/skills/windsurf-load-scale/SKILL.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/windsurf-pack/skills/windsurf-load-scale/SKILL.md` |
| **Category** | development |
| **Subcategory** | tools |
| **Tags** | development |
| **Created** | 2026-01-06 |
| **Updated** | 2026-01-06 |
| **File Hash** | `232935490e9aee29...` |

## Description

export default function () {
  const response = http.post(
    'https://api.windsurf.com/v1/resource',
    JSON.stringify({ test: true }),
    {
      headers: {
        'ContentType': 'application/json',
        'Authorization': Bearer ${__ENV.WINDSURF_API_KEY},
      },
    }
  );

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/windsurf-pack/skills/windsurf-load-scale/SKILL.md)*
