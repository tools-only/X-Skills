# Drm

| Property | Value |
|----------|-------|
| **Name** | Drm |
| **Repository** | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/drm.md) (⭐ 612) |
| **Original Path** | `kernel/subsystem/drm.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-12 |
| **Updated** | 2026-02-12 |
| **File Hash** | `c07f33916e2b0818...` |

## Description

Calling sleeping functions from atomic context causes kernel warnings, system
instability, and potential deadlocks. DRM/KMS display drivers have multiple
code paths that execute in atomic context where sleeping is forbidden.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/drm.md)*
