# Mm Vma

| Property | Value |
|----------|-------|
| **Name** | Mm Vma |
| **Repository** | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/mm-vma.md) (⭐ 629) |
| **Original Path** | `kernel/subsystem/mm-vma.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-13 |
| **Updated** | 2026-02-16 |
| **File Hash** | `efe0703128a10fe2...` |

## Description

Dereferencing a parent/owner pointer from a SLAB_TYPESAFE_BY_RCU object
after dropping the object's refcount causes useafterfree when the object has
been recycled to a different owner. The owner can exit and free its backing
structure in the window between the refcount drop and the dereference.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/mm-vma.md)*
