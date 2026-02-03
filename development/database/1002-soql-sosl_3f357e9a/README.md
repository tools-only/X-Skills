# SOQL and SOSL

| Property | Value |
|----------|-------|
| **Name** | SOQL and SOSL |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/salesforce-developer/references/soql-sosl.md) (⭐ 216) |
| **Original Path** | `skills/salesforce-developer/references/soql-sosl.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-01-19 |
| **Updated** | 2026-01-29 |
| **File Hash** | `3f357e9a104cd492...` |

## Description

sql
SELECT Id, Name, Industry, AnnualRevenue
FROM Account
WHERE Industry = 'Technology'
AND AnnualRevenue > 1000000
ORDER BY AnnualRevenue DESC NULLS LAST
LIMIT 100
OFFSET 0

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/salesforce-developer/references/soql-sosl.md)*
