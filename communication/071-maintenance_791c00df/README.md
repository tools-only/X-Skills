# Database Maintenance

| Property | Value |
|----------|-------|
| **Name** | Database Maintenance |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/postgres-pro/references/maintenance.md) (⭐ 216) |
| **Original Path** | `skills/postgres-pro/references/maintenance.md` |
| **Category** | communication |
| **Subcategory** | email |
| **Tags** | communication |
| **Created** | 2025-12-15 |
| **Updated** | 2026-01-29 |
| **File Hash** | `791c00df8abc92f7...` |

## Description

PostgreSQL uses MVCC (MultiVersion Concurrency Control):
 Updates/deletes don't remove old rows immediately
 Old rows marked as "dead tuples"
 VACUUM reclaims space from dead tuples
 Without VACUUM: table bloat, degraded performance, transaction ID wraparound

**Tags:** `communication`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/postgres-pro/references/maintenance.md)*
