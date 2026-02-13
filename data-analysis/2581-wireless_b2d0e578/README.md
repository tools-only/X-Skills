# Wireless

| Property | Value |
|----------|-------|
| **Name** | Wireless |
| **Repository** | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/wireless.md) (⭐ 612) |
| **Original Path** | `kernel/subsystem/wireless.md` |
| **Category** | data-analysis |
| **Subcategory** | processing |
| **Tags** | data analysis |
| **Created** | 2026-02-12 |
| **Updated** | 2026-02-12 |
| **File Hash** | `b2d0e5788966ba26...` |

## Description

Incorrect placement of BSS_CHANGED_ event handling between mac80211
callbacks causes WARN_ON_ONCE splats at runtime, functional regressions
(association failures, beacon processing errors), or silent hardware
misbehavior that is difficult to diagnose.

**Tags:** `data analysis`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/wireless.md)*
