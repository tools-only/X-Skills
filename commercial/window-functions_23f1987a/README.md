# Window Functions

| Property | Value |
|----------|-------|
| **Name** | Window Functions |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/sql-pro/references/window-functions.md) (⭐ 216) |
| **Original Path** | `skills/sql-pro/references/window-functions.md` |
| **Category** | commercial |
| **Subcategory** | ecommerce |
| **Tags** | commercial |
| **Created** | 2025-12-15 |
| **Updated** | 2026-01-29 |
| **File Hash** | `23f1987aa8223d25...` |

## Description

sql
 ROW_NUMBER: Sequential numbering within partition
SELECT
    customer_id,
    order_date,
    total,
    ROW_NUMBER() OVER (PARTITION BY customer_id ORDER BY order_date DESC) as row_num
FROM orders;

**Tags:** `commercial`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/sql-pro/references/window-functions.md)*
