# Performance Optimization

| Property | Value |
|----------|-------|
| **Name** | Performance Optimization |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/postgres-pro/references/performance.md) (⭐ 216) |
| **Original Path** | `skills/postgres-pro/references/performance.md` |
| **Category** | content-creation |
| **Subcategory** | writing |
| **Tags** | content creation |
| **Created** | 2025-12-15 |
| **Updated** | 2026-01-29 |
| **File Hash** | `61f5b3c5370c59f6...` |

## Description

sql
 Basic EXPLAIN ANALYZE
EXPLAIN (ANALYZE, BUFFERS, VERBOSE)
SELECT u.id, u.name, COUNT(o.id) as order_count
FROM users u
LEFT JOIN orders o ON u.id = o.user_id
WHERE u.created_at > '20240101'
GROUP BY u.id, u.name;

**Tags:** `content creation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/postgres-pro/references/performance.md)*
