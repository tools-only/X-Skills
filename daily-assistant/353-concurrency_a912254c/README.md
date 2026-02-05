# Concurrency

| Property | Value |
|----------|-------|
| **Name** | Concurrency |
| **Repository** | [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/skills/swiftdata-code-review/references/concurrency.md) (⭐ 17) |
| **Original Path** | `skills/swiftdata-code-review/references/concurrency.md` |
| **Category** | daily-assistant |
| **Subcategory** | tasks |
| **Tags** | daily assistant |
| **Created** | 2026-01-11 |
| **Updated** | 2026-01-11 |
| **File Hash** | `a912254c650debbb...` |

## Description

swift
@ModelActor
actor DataHandler {
    func importItems(_ data: [ImportData]) throws {
        for item in data {
            modelContext.insert(Item(name: item.name))
        }
        try modelContext.save()
    }

**Tags:** `daily assistant`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/skills/swiftdata-code-review/references/concurrency.md)*
