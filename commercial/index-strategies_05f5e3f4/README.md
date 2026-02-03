# Index Strategies

| Property | Value |
|----------|-------|
| **Name** | Index Strategies |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/database-optimizer/references/index-strategies.md) (⭐ 216) |
| **Original Path** | `skills/database-optimizer/references/index-strategies.md` |
| **Category** | commercial |
| **Subcategory** | ecommerce |
| **Tags** | commercial |
| **Created** | 2025-12-15 |
| **Updated** | 2026-01-29 |
| **File Hash** | `05f5e3f4e603ca69...` |

## Description

sql
 PostgreSQL: Find queries missing indexes
SELECT query, calls, total_exec_time, mean_exec_time
FROM pg_stat_statements
WHERE mean_exec_time > 100
ORDER BY total_exec_time DESC
LIMIT 20;

**Tags:** `commercial`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/database-optimizer/references/index-strategies.md)*
