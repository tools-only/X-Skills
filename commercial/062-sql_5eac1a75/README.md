# Dialect Differences

| Property | Value |
|----------|-------|
| **Name** | Dialect Differences |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/sql-pro/references/dialect-differences.md) (⭐ 216) |
| **Original Path** | `skills/sql-pro/references/dialect-differences.md` |
| **Category** | commercial |
| **Subcategory** | ecommerce |
| **Tags** | commercial |
| **Created** | 2025-12-15 |
| **Updated** | 2026-01-29 |
| **File Hash** | `5eac1a75d44e65da...` |

## Description

sql
 PostgreSQL
CREATE TABLE users (
    user_id SERIAL PRIMARY KEY,   or BIGSERIAL for BIGINT
    name VARCHAR(100)
);
 Alternative (PostgreSQL 10+)
CREATE TABLE users (
    user_id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    name VARCHAR(100)
);

**Tags:** `commercial`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/sql-pro/references/dialect-differences.md)*
