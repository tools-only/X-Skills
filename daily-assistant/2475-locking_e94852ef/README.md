# Locking

| Property | Value |
|----------|-------|
| **Name** | Locking |
| **Repository** | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/locking.md) (⭐ 612) |
| **Original Path** | `kernel/subsystem/locking.md` |
| **Category** | daily-assistant |
| **Subcategory** | notes |
| **Tags** | daily assistant |
| **Created** | 2026-02-11 |
| **Updated** | 2026-02-12 |
| **File Hash** | `e94852efca64d11d...` |

## Description

Using the wrong lock type for the execution context causes deadlocks (sleeping
in atomic context), missed wakeups, or priority inversion. The table below
shows which lock types provide correct protection for data shared across the
listed contexts.

**Tags:** `daily assistant`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/locking.md)*
