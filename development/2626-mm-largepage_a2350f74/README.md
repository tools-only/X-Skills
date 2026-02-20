# Mm Largepage

| Property | Value |
|----------|-------|
| **Name** | Mm Largepage |
| **Repository** | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/mm-largepage.md) (⭐ 629) |
| **Original Path** | `kernel/subsystem/mm-largepage.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-13 |
| **Updated** | 2026-02-13 |
| **File Hash** | `a2350f74bf0d4cf4...` |

## Description

Misunderstanding which flags are perpage vs perfolio leads to bugs where code
checks or sets state on the wrong struct page. A common mistake is assuming all
flags work like small pages when operating on subpages of a large folio.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/mm-largepage.md)*
