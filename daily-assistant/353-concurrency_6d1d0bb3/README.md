# Concurrency

| Property | Value |
|----------|-------|
| **Name** | Concurrency |
| **Repository** | [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/skills/swift-code-review/references/concurrency.md) (⭐ 17) |
| **Original Path** | `skills/swift-code-review/references/concurrency.md` |
| **Category** | daily-assistant |
| **Subcategory** | tasks |
| **Tags** | daily assistant |
| **Created** | 2026-01-11 |
| **Updated** | 2026-01-11 |
| **File Hash** | `6d1d0bb332db12c5...` |

## Description

swift
// BAD  sequential awaits on independent operations
let user = await fetchUser()
let avatar = await fetchAvatar(for: user.id)      // Waits unnecessarily
let prefs = await fetchPreferences(for: user.id)  // Waits unnecessarily

**Tags:** `daily assistant`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/skills/swift-code-review/references/concurrency.md)*
