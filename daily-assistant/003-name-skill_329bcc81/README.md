# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/skills/elixir-performance-review/SKILL.md) (⭐ 17) |
| **Original Path** | `skills/elixir-performance-review/SKILL.md` |
| **Category** | daily-assistant |
| **Subcategory** | tasks |
| **Tags** | daily assistant |
| **Created** | 2026-02-05 |
| **Updated** | 2026-02-05 |
| **File Hash** | `329bcc81cc5d5a62...` |

## Description

Single GenServer for lowthroughput  Not all state needs horizontal scaling
 Synchronous calls for critical paths  Consistency may require it
 Inmemory state without ETS  ETS has overhead for small state
 Enum over Stream for small collections  Stream overhead not worth it

**Tags:** `daily assistant`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/skills/elixir-performance-review/SKILL.md)*
