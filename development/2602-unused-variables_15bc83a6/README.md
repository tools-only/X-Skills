# Unused Variables

| Property | Value |
|----------|-------|
| **Name** | Unused Variables |
| **Repository** | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/unused-variables.md) (⭐ 90) |
| **Original Path** | `plugins/scv-scan/skills/scv-scan/references/unused-variables.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-13 |
| **Updated** | 2026-02-13 |
| **File Hash** | `15bc83a6b69e83bd...` |

## Description

function deposit(uint256 amount, bytes memory data) external {
        // data parameter never used — possible missing validation
        totalDeposits += amount;
    }

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/unused-variables.md)*
