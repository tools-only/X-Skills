# Limitations of Dynamic SQL

| Property | Value |
|----------|-------|
| **Name** | Limitations of Dynamic SQL |
| **Repository** | [dbt-labs/dbt-agent-skills](https://raw.githubusercontent.com/dbt-labs/dbt-agent-skills/main/skills/migrating-dbt-core-to-fusion/manual_fixes/dynamic_sql.md) (⭐ 10) |
| **Original Path** | `skills/migrating-dbt-core-to-fusion/manual_fixes/dynamic_sql.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-01-29 |
| **Updated** | 2026-01-29 |
| **File Hash** | `94831248c174055c...` |

## Description

Some SQL patterns that work in legacy dbt may not be supported by Fusion's static analysis, such as:
sql
PIVOT(...) FOR column_name IN (ANY)

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [dbt-labs/dbt-agent-skills](https://raw.githubusercontent.com/dbt-labs/dbt-agent-skills/main/skills/migrating-dbt-core-to-fusion/manual_fixes/dynamic_sql.md)*
