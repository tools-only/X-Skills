# Networking

| Property | Value |
|----------|-------|
| **Name** | Networking |
| **Repository** | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/networking.md) (⭐ 629) |
| **Original Path** | `kernel/subsystem/networking.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-11 |
| **Updated** | 2026-02-14 |
| **File Hash** | `44c7f9d8c5aa3a5c...` |

## Description

skb_put(), skb_push(), and skb_pull() modify the data boundaries of a
socket buffer. Passing untrusted or unchecked lengths causes a kernel panic
(DoS). The bounds checks fire before memory is corrupted, so the result is a
crash rather than a silent overflow, but it is still a bug.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/networking.md)*
