# Caveats for Redshift

| Property | Value |
|----------|-------|
| **Name** | Caveats for Redshift |
| **Repository** | [dbt-labs/dbt-agent-skills](https://raw.githubusercontent.com/dbt-labs/dbt-agent-skills/main/skills/adding-dbt-unit-test/warehouses/redshift/caveats.md) (⭐ 10) |
| **Original Path** | `skills/adding-dbt-unit-test/warehouses/redshift/caveats.md` |
| **Category** | development |
| **Subcategory** | testing |
| **Tags** | development |
| **Created** | 2026-01-29 |
| **Updated** | 2026-01-29 |
| **File Hash** | `5bec6b91712be0b9...` |

## Description

Redshift doesn't support unit tests when the SQL in the common table expression (CTE) contains functions such as LISTAGG, MEDIAN, PERCENTILE_CONT, and so on. These functions must be executed against a usercreated table. dbt combines given rows to be part of the CTE, which Redshift does not support.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [dbt-labs/dbt-agent-skills](https://raw.githubusercontent.com/dbt-labs/dbt-agent-skills/main/skills/adding-dbt-unit-test/warehouses/redshift/caveats.md)*
