# Setup

| Property | Value |
|----------|-------|
| **Name** | Setup |
| **Repository** | [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/skills/sqlite-vec/references/setup.md) (⭐ 17) |
| **Original Path** | `skills/sqlite-vec/references/setup.md` |
| **Category** | development |
| **Subcategory** | tools |
| **Tags** | development |
| **Created** | 2025-12-21 |
| **Updated** | 2025-12-21 |
| **File Hash** | `2a6ed82347c51ab8...` |

## Description

db = sqlite3.connect(":memory:")
db.enable_load_extension(True)
sqlite_vec.load(db)
db.enable_load_extension(False)

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/skills/sqlite-vec/references/setup.md)*
