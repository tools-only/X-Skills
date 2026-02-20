# Mm Alloc

| Property | Value |
|----------|-------|
| **Name** | Mm Alloc |
| **Repository** | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/mm-alloc.md) (⭐ 629) |
| **Original Path** | `kernel/subsystem/mm-alloc.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-13 |
| **Updated** | 2026-02-14 |
| **File Hash** | `93875e44bf617608...` |

## Description

Using the wrong GFP flag causes sleeping in atomic context (deadlock/BUG),
filesystem or IO recursion (deadlock), or silent allocation failures when the
caller assumes success. Verify the allocation context matches the flag.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/mm-alloc.md)*
