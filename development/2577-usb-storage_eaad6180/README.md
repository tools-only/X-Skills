# Usb Storage

| Property | Value |
|----------|-------|
| **Name** | Usb Storage |
| **Repository** | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/usb-storage.md) (⭐ 612) |
| **Original Path** | `kernel/subsystem/usb-storage.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-12 |
| **Updated** | 2026-02-12 |
| **File Hash** | `eaad618009d0d8c8...` |

## Description

Specifying unnecessary subclass or protocol overrides in UNUSUAL_DEV() entries
causes the driver to emit a dev_notice on every device insertion, asking users
to report the unneeded entry (see get_device_info() in
drivers/usb/storage/usb.c). This creates persistent log noise for end users.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/usb-storage.md)*
