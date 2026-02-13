# Kho

| Property | Value |
|----------|-------|
| **Name** | Kho |
| **Repository** | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/kho.md) (⭐ 612) |
| **Original Path** | `kernel/subsystem/kho.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-12 |
| **Updated** | 2026-02-12 |
| **File Hash** | `6a325a13ed1c631e...` |

## Description

c
// WRONG: Missing enabled check
static int __init my_kho_init(void)
{
    err = kho_add_subtree("my_node", fdt);
    // NULL deref on kho_out.fdt if KHO is disabled
}

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/kho.md)*
