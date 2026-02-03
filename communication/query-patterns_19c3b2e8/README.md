# Query Patterns

| Property | Value |
|----------|-------|
| **Name** | Query Patterns |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/sql-pro/references/query-patterns.md) (⭐ 216) |
| **Original Path** | `skills/sql-pro/references/query-patterns.md` |
| **Category** | communication |
| **Subcategory** | email |
| **Tags** | communication |
| **Created** | 2025-12-15 |
| **Updated** | 2026-01-29 |
| **File Hash** | `19c3b2e803336b5d...` |

## Description

Bill of materials (parts explosion)
WITH RECURSIVE parts_explosion AS (
    SELECT
        part_id,
        component_id,
        quantity,
        1 as level,
        ARRAY[part_id] as path
    FROM bill_of_materials
    WHERE part_id = 'PRODUCT123'

**Tags:** `communication`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/sql-pro/references/query-patterns.md)*
