# Caching Strategy

| Property | Value |
|----------|-------|
| **Name** | Caching Strategy |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/supabase-pack/skills/supabase-performance-tuning/references/caching-strategy.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/supabase-pack/skills/supabase-performance-tuning/references/caching-strategy.md` |
| **Category** | data-analysis |
| **Subcategory** | processing |
| **Tags** | data analysis |
| **Created** | 2026-01-06 |
| **Updated** | 2026-01-06 |
| **File Hash** | `e7386a4845ba5be1...` |

## Description

const cache = new LRUCache<string, any>({
  max: 1000,
  ttl: 30000, // 1 minute
  updateAgeOnGet: true,
});

**Tags:** `data analysis`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/supabase-pack/skills/supabase-performance-tuning/references/caching-strategy.md)*
