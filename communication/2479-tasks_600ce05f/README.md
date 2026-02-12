# Tasks

| Property | Value |
|----------|-------|
| **Name** | Tasks |
| **Repository** | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-user-scoped-features/tasks.md) (⭐ 19) |
| **Original Path** | `openspec/changes/add-user-scoped-features/tasks.md` |
| **Category** | communication |
| **Subcategory** | messaging |
| **Tags** | communication |
| **Created** | 2026-01-21 |
| **Updated** | 2026-01-21 |
| **File Hash** | `600ce05f93c07ae5...` |

## Description

[x] 1.1 创建 chat_sessions 表 (session_id, user_id, title, created_at, updated_at, last_message_at, message_count)
 [x] 1.2 创建 chat_messages 表 (id, session_id, user_id, role, content, metadata, created_at)
 [x] 1.3 users 表增加 is_admin 字段 (UInt8, 默认0)
 [x] 1.4 更新 sync_history 表增加 user_id 字段

**Tags:** `communication`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-user-scoped-features/tasks.md)*
