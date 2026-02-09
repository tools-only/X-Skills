# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/clerk-pack/skills/clerk-debug-bundle/SKILL.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/clerk-pack/skills/clerk-debug-bundle/SKILL.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-01-09 |
| **Updated** | 2026-01-09 |
| **File Hash** | `390692f83abce2e7...` |

## Description

// Get package versions
  try {
    const pkg = require('../package.json')
    debug.packages = {
      '@clerk/nextjs': pkg.dependencies?.['@clerk/nextjs'],
      '@clerk/clerkreact': pkg.dependencies?.['@clerk/clerkreact'],
      'next': pkg.dependencies?.['next']
    }
  } catch {}

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/clerk-pack/skills/clerk-debug-bundle/SKILL.md)*
