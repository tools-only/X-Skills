# Mysql Stored Procedure Best Practices

| Property | Value |
|----------|-------|
| **Name** | Mysql Stored Procedure Best Practices |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/database/stored-procedure-generator/skills/generating-stored-procedures/references/mysql_stored_procedure_best_practices.md) (⭐ 1.3k) |
| **Original Path** | `plugins/database/stored-procedure-generator/skills/generating-stored-procedures/references/mysql_stored_procedure_best_practices.md` |
| **Category** | commercial |
| **Subcategory** | ecommerce |
| **Tags** | commercial |
| **Created** | 2026-01-20 |
| **Updated** | 2026-01-20 |
| **File Hash** | `3ab757cf41a5cd91...` |

## Description

sql
DELIMITER //
CREATE PROCEDURE procedure_name(
    IN param1 INT,
    OUT param2 VARCHAR(255),
    INOUT param3 DECIMAL(10,2)
)
BEGIN
     Procedure body
END //
DELIMITER ;

**Tags:** `commercial`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/database/stored-procedure-generator/skills/generating-stored-procedures/references/mysql_stored_procedure_best_practices.md)*
