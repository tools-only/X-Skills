# Mm Folio

| Property | Value |
|----------|-------|
| **Name** | Mm Folio |
| **Repository** | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/mm-folio.md) (⭐ 629) |
| **Original Path** | `kernel/subsystem/mm-folio.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-13 |
| **Updated** | 2026-02-13 |
| **File Hash** | `258ed4d2660b981a...` |

## Description

Lazyfree folios (!folio_test_swapbacked()) may be discarded without writeback
if clean. The reclaim path must check dirty status before refcount, and only
call folio_set_swapbacked() when the folio is genuinely dirty:

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/mm-folio.md)*
