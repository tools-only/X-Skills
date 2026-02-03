# .eventOccurred,

| Property | Value |
|----------|-------|
| **Name** | .eventOccurred, |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/swift-expert/references/async-concurrency.md) (⭐ 216) |
| **Original Path** | `skills/swift-expert/references/async-concurrency.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2025-12-15 |
| **Updated** | 2026-01-29 |
| **File Hash** | `7555824cab883df0...` |

## Description

swift
// Async function
func fetchUser(id: Int) async throws > User {
    let url = URL(string: "https://api.example.com/users/\(id)")!
    let (data, _) = try await URLSession.shared.data(from: url)
    return try JSONDecoder().decode(User.self, from: data)
}

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/swift-expert/references/async-concurrency.md)*
