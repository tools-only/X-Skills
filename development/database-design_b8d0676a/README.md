# Database Design

| Property | Value |
|----------|-------|
| **Name** | Database Design |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/sql-pro/references/database-design.md) (⭐ 216) |
| **Original Path** | `skills/sql-pro/references/database-design.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2025-12-15 |
| **Updated** | 2026-01-29 |
| **File Hash** | `b8d0676ac9121399...` |

## Description

sql
 1NF: Atomic values, no repeating groups
 Bad: Nonatomic phone column
CREATE TABLE customers_bad (
    customer_id INT PRIMARY KEY,
    name VARCHAR(100),
    phones VARCHAR(500)   "5551234,5555678,5559012"
);

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/sql-pro/references/database-design.md)*
