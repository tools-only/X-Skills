# Networking

| Property | Value |
|----------|-------|
| **Name** | Networking |
| **Repository** | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/networking.md) (⭐ 611) |
| **Original Path** | `kernel/subsystem/networking.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-11 |
| **Updated** | 2026-02-11 |
| **File Hash** | `a82d52ca143e015a...` |

## Description

skb_put(), skb_push(), and skb_pull() modify the data boundaries of a
socket buffer. Each operation has a bounds check that triggers a kernel panic
on violation:

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/networking.md)*
