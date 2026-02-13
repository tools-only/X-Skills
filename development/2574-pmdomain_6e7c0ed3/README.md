# Pmdomain

| Property | Value |
|----------|-------|
| **Name** | Pmdomain |
| **Repository** | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/pmdomain.md) (⭐ 612) |
| **Original Path** | `kernel/subsystem/pmdomain.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-12 |
| **Updated** | 2026-02-12 |
| **File Hash** | `6e7c0ed3fb8d0eec...` |

## Description

The stay_on flag is cleared when the provider's sync_state callback runs.
For OFbased providers, this is handled by genpd_provider_sync_state() or
of_genpd_sync_state(), both in drivers/pmdomain/core.c. These functions
set genpd>stay_on = false and then attempt genpd_power_off().

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/pmdomain.md)*
