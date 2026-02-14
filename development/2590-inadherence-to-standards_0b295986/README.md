# Inadherence To Standards

| Property | Value |
|----------|-------|
| **Name** | Inadherence To Standards |
| **Repository** | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/inadherence-to-standards.md) (⭐ 90) |
| **Original Path** | `plugins/scv-scan/skills/scv-scan/references/inadherence-to-standards.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-13 |
| **Updated** | 2026-02-13 |
| **File Hash** | `0b295986e5907c83...` |

## Description

function depositToken(IERC20 token, uint256 amount) external {
    uint256 balBefore = token.balanceOf(address(this));
    token.safeTransferFrom(msg.sender, address(this), amount);
    uint256 received = token.balanceOf(address(this))  balBefore;
    deposits[msg.sender] += received;
}

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/inadherence-to-standards.md)*
