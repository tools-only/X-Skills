# Wait for specific standby

| Property | Value |
|----------|-------|
| **Name** | Wait for specific standby |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/postgres-pro/references/replication.md) (⭐ 216) |
| **Original Path** | `skills/postgres-pro/references/replication.md` |
| **Category** | commercial |
| **Subcategory** | ecommerce |
| **Tags** | commercial |
| **Created** | 2025-12-15 |
| **Updated** | 2026-01-29 |
| **File Hash** | `c4993a7c459f10f4...` |

## Description

sql
 postgresql.conf
wal_level = replica
max_wal_senders = 10
max_replication_slots = 10
wal_keep_size = 1GB   Or 1024MB for older versions
hot_standby = on
archive_mode = on
archive_command = 'cp %p /var/lib/postgresql/wal_archive/%f'

**Tags:** `commercial`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/postgres-pro/references/replication.md)*
