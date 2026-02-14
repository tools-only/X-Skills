# Unchecked Return Values

| Property | Value |
|----------|-------|
| **Name** | Unchecked Return Values |
| **Repository** | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/unchecked-return-values.md) (⭐ 90) |
| **Original Path** | `plugins/scv-scan/skills/scv-scan/references/unchecked-return-values.md` |
| **Category** | commercial |
| **Subcategory** | business |
| **Tags** | commercial |
| **Created** | 2026-02-13 |
| **Updated** | 2026-02-13 |
| **File Hash** | `ae8ce7f878240113...` |

## Description

function payout(address to, uint256 amount) external {
    // .call() return value captured but never checked
    (bool success,) = to.call{value: amount}("");
    // success could be false, but execution continues
    totalPaid += amount;
}

**Tags:** `commercial`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/unchecked-return-values.md)*
