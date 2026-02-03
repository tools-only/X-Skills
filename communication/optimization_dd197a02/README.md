# Query Optimization

| Property | Value |
|----------|-------|
| **Name** | Query Optimization |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/sql-pro/references/optimization.md) (⭐ 216) |
| **Original Path** | `skills/sql-pro/references/optimization.md` |
| **Category** | communication |
| **Subcategory** | email |
| **Tags** | communication |
| **Created** | 2025-12-15 |
| **Updated** | 2026-01-29 |
| **File Hash** | `dd197a02491cda9e...` |

## Description

MySQL EXPLAIN
EXPLAIN FORMAT=JSON
SELECT  FROM orders o
INNER JOIN customers c ON o.customer_id = c.customer_id
WHERE o.order_date >= '20240101'
  AND c.country = 'US';

**Tags:** `communication`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/sql-pro/references/optimization.md)*
