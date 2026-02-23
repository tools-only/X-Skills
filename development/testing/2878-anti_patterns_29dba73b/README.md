# Anti Patterns

| Property | Value |
|----------|-------|
| **Name** | Anti Patterns |
| **Repository** | [glebis/claude-skills](https://raw.githubusercontent.com/glebis/claude-skills/main/tdd/references/anti_patterns.md) (⭐ 18) |
| **Original Path** | `tdd/references/anti_patterns.md` |
| **Category** | development |
| **Subcategory** | testing |
| **Tags** | development |
| **Created** | 2026-02-19 |
| **Updated** | 2026-02-19 |
| **File Hash** | `29dba73b06e0731f...` |

## Description

Fix: Test through public interfaces only. Assert on return values, side effects visible to callers, or observable state changes.
Why it matters: Implementationdetail tests break on every refactor, even when behavior is preserved. They test HOW the code works, not WHAT it does.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [glebis/claude-skills](https://raw.githubusercontent.com/glebis/claude-skills/main/tdd/references/anti_patterns.md)*
