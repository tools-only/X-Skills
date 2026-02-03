# Migration Gotchas & Common Pitfalls

| Property | Value |
|----------|-------|
| **Name** | Migration Gotchas & Common Pitfalls |
| **Repository** | [adityamparikh/maven-to-gradle](https://raw.githubusercontent.com/adityamparikh/maven-to-gradle/main/references/gotchas.md) (⭐ 10) |
| **Original Path** | `references/gotchas.md` |
| **Category** | development |
| **Subcategory** | testing |
| **Tags** | development |
| **Created** | 2026-02-02 |
| **Updated** | 2026-02-02 |
| **File Hash** | `66e17050ddf98e62...` |

## Description

Critical difference: Maven's compile scope is transitive. Gradle's implementation is NOT. If downstream modules need a dependency transitively, use the javalibrary plugin and api configuration.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [adityamparikh/maven-to-gradle](https://raw.githubusercontent.com/adityamparikh/maven-to-gradle/main/references/gotchas.md)*
