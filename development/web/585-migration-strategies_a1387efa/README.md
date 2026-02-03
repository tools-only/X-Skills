# str

| Property | Value |
|----------|-------|
| **Name** | str |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/legacy-modernizer/references/migration-strategies.md) (⭐ 216) |
| **Original Path** | `skills/legacy-modernizer/references/migration-strategies.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2025-12-15 |
| **Updated** | 2026-01-29 |
| **File Hash** | `a1387efa79f4970a...` |

## Description

python
 Phase 1: Dual write to both databases
class DualWriteUserRepository:
    def __init__(self, legacy_db, modern_db: AsyncSession):
        self.legacy = legacy_db
        self.modern = modern_db

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/legacy-modernizer/references/migration-strategies.md)*
