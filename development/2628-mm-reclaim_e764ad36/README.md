# Mm Reclaim

| Property | Value |
|----------|-------|
| **Name** | Mm Reclaim |
| **Repository** | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/mm-reclaim.md) (⭐ 629) |
| **Original Path** | `kernel/subsystem/mm-reclaim.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-13 |
| **Updated** | 2026-02-14 |
| **File Hash** | `e764ad369f9f0a60...` |

## Description

Incorrect tag handling causes data loss (dirty pages skipped during sync) or
writeback livelock (sync never completes because new dirty pages keep appearing).
Review any code that starts writeback or implements >writepages.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/mm-reclaim.md)*
