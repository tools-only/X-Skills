# Timestamp Dependence

| Property | Value |
|----------|-------|
| **Name** | Timestamp Dependence |
| **Repository** | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/timestamp-dependence.md) (⭐ 90) |
| **Original Path** | `plugins/scv-scan/skills/scv-scan/references/timestamp-dependence.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-13 |
| **Updated** | 2026-02-13 |
| **File Hash** | `0f401084feb4256d...` |

## Description

// Tight time window vulnerable to manipulation
function claimBonus() external {
    // 15second window — validator can push timestamp to include/exclude
    require(block.timestamp >= deadline && block.timestamp <= deadline + 15);
    _sendBonus(msg.sender);
}

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/timestamp-dependence.md)*
