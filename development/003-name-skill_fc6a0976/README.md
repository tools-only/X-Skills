# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/lokalise-pack/skills/lokalise-rate-limits/SKILL.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/lokalise-pack/skills/lokalise-rate-limits/SKILL.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-01-26 |
| **Updated** | 2026-01-26 |
| **File Hash** | `fc6a0976f55975e8...` |

## Description

// Lokalise: 6 req/sec, 10 concurrent per project
const queue = new PQueue({
  concurrency: 5,        // Conservative concurrent limit
  interval: 1000,        // 1 second
  intervalCap: 5,        // Max 5 per second (leave headroom)
  carryoverConcurrencyCount: true,
});

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/lokalise-pack/skills/lokalise-rate-limits/SKILL.md)*
