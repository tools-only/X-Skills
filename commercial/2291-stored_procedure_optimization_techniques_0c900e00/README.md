# Stored Procedure Optimization Techniques

| Property | Value |
|----------|-------|
| **Name** | Stored Procedure Optimization Techniques |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/database/stored-procedure-generator/skills/generating-stored-procedures/references/stored_procedure_optimization_techniques.md) (⭐ 1.3k) |
| **Original Path** | `plugins/database/stored-procedure-generator/skills/generating-stored-procedures/references/stored_procedure_optimization_techniques.md` |
| **Category** | commercial |
| **Subcategory** | ecommerce |
| **Tags** | commercial |
| **Created** | 2026-01-20 |
| **Updated** | 2026-01-20 |
| **File Hash** | `0c900e00eabcd937...` |

## Description

sql
 SLOW: Rowbyrow cursor processing
CREATE FUNCTION update_prices_slow()
RETURNS VOID AS $$
DECLARE
    rec RECORD;
BEGIN
    FOR rec IN SELECT id, price FROM products LOOP
        UPDATE products SET price = rec.price  1.1 WHERE id = rec.id;
    END LOOP;
END;
$$ LANGUAGE plpgsql;

**Tags:** `commercial`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/database/stored-procedure-generator/skills/generating-stored-procedures/references/stored_procedure_optimization_techniques.md)*
