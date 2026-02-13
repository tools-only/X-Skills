# Power Management Subsystem Details

| Property | Value |
|----------|-------|
| **Name** | Power Management Subsystem Details |
| **Repository** | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/pm.md) (⭐ 612) |
| **Original Path** | `kernel/subsystem/pm.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-12 |
| **Updated** | 2026-02-12 |
| **File Hash** | `5127098d16f1a0f7...` |

## Description

Misunderstanding which return values are possible from each Runtime PM API
leads to incorrect error handling, missed wakeups, or treating success as
failure. The three base functions have distinct return semantics that callers
must respect.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/pm.md)*
