# Memory Management Subsystem Details

| Property | Value |
|----------|-------|
| **Name** | Memory Management Subsystem Details |
| **Repository** | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/mm.md) (⭐ 612) |
| **Original Path** | `kernel/subsystem/mm.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-11 |
| **Updated** | 2026-02-12 |
| **File Hash** | `409e957e73e3c2da...` |

## Description

Incorrect PTE flag combinations cause data corruption (dirty data silently
dropped), security holes (writable pages that should be readonly), and kernel
crashes on architectures that trap invalid combinations. Review any code that
constructs or modifies PTEs for these invariants.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/mm.md)*
