# Scheduler

| Property | Value |
|----------|-------|
| **Name** | Scheduler |
| **Repository** | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/scheduler.md) (⭐ 611) |
| **Original Path** | `kernel/subsystem/scheduler.md` |
| **Category** | daily-assistant |
| **Subcategory** | scheduling |
| **Tags** | daily assistant |
| **Created** | 2026-02-11 |
| **Updated** | 2026-02-11 |
| **File Hash** | `f3365e4e6a24ac45...` |

## Description

set_current_state() includes a barrier so the state write is ordered
  relative to subsequent memory accesses (the condition check)
 Voluntary sleep pattern: set state BEFORE checking the condition, otherwise
  a wakeup between the check and the state change is lost

**Tags:** `daily assistant`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/scheduler.md)*
