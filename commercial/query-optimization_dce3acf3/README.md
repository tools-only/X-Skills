# Query Optimization

| Property | Value |
|----------|-------|
| **Name** | Query Optimization |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/database-optimizer/references/query-optimization.md) (⭐ 216) |
| **Original Path** | `skills/database-optimizer/references/query-optimization.md` |
| **Category** | commercial |
| **Subcategory** | ecommerce |
| **Tags** | commercial |
| **Created** | 2025-12-15 |
| **Updated** | 2026-01-29 |
| **File Hash** | `dce3acf3b12dee59...` |

## Description

sql
 Get actual execution statistics
EXPLAIN (ANALYZE, BUFFERS, VERBOSE, TIMING)
SELECT u.id, u.name, COUNT(o.id) as order_count
FROM users u
LEFT JOIN orders o ON u.id = o.user_id
WHERE u.created_at > NOW()  INTERVAL '30 days'
GROUP BY u.id, u.name
HAVING COUNT(o.id) > 5;

**Tags:** `commercial`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/database-optimizer/references/query-optimization.md)*
