# Requirement Violation

| Property | Value |
|----------|-------|
| **Name** | Requirement Violation |
| **Repository** | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/requirement-violation.md) (⭐ 90) |
| **Original Path** | `plugins/scv-scan/skills/scv-scan/references/requirement-violation.md` |
| **Category** | commercial |
| **Subcategory** | business |
| **Tags** | commercial |
| **Created** | 2026-02-13 |
| **Updated** | 2026-02-13 |
| **File Hash** | `8c2a25b18d7ed594...` |

## Description

// External contract assumption mismatch
function processPayment(IERC20 token, uint256 amount) external {
    // Assumes transfer returns true, but some tokens don't return a value
    // require reverts on tokens like USDT
    require(token.transfer(msg.sender, amount), "failed");
}

**Tags:** `commercial`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/requirement-violation.md)*
