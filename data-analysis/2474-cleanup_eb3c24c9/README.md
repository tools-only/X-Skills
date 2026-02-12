# Cleanup

| Property | Value |
|----------|-------|
| **Name** | Cleanup |
| **Repository** | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/cleanup.md) (⭐ 611) |
| **Original Path** | `kernel/subsystem/cleanup.md` |
| **Category** | data-analysis |
| **Subcategory** | processing |
| **Tags** | data analysis |
| **Created** | 2026-02-11 |
| **Updated** | 2026-02-11 |
| **File Hash** | `eb3c24c90abf7b8c...` |

## Description

Mitigation patterns:
 Setting the variable to NULL before early return
 Using no_free_ptr() before early return to inhibit cleanup
 Using a cleanup function that checks IS_ERR_OR_NULL()

**Tags:** `data analysis`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/cleanup.md)*
