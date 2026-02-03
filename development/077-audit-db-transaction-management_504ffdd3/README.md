# Audit Db Transaction Management

| Property | Value |
|----------|-------|
| **Name** | Audit Db Transaction Management |
| **Repository** | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/llms/audit-db-transaction-management.md) (⭐ 3.2k) |
| **Original Path** | `llms/audit-db-transaction-management.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-03 |
| **Updated** | 2026-02-03 |
| **File Hash** | `504ffdd3fc5257ce...` |

## Description

This prompt guides you through auditing the MCP Gateway codebase for proper database transaction management. The goal is to identify endpoints missing explicit db.commit(); db.close() calls, which cause PostgreSQL connection leaks under load.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/llms/audit-db-transaction-management.md)*
