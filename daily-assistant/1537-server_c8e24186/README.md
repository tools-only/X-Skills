# Server

| Property | Value |
|----------|-------|
| **Name** | Server |
| **Repository** | [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/skills/wish-ssh-code-review/references/server.md) (⭐ 17) |
| **Original Path** | `skills/wish-ssh-code-review/references/server.md` |
| **Category** | daily-assistant |
| **Subcategory** | notes |
| **Tags** | daily assistant |
| **Created** | 2025-12-21 |
| **Updated** | 2025-12-21 |
| **File Hash** | `c8e24186db4ac55e...` |

## Description

go
// BAD  generates new key each start (fingerprint changes)
s, err := wish.NewServer(
    wish.WithAddress(":22"),
    // no host key specified  generates random
)

**Tags:** `daily assistant`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/skills/wish-ssh-code-review/references/server.md)*
