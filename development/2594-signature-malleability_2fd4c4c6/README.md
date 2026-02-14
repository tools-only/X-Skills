# Signature Malleability

| Property | Value |
|----------|-------|
| **Name** | Signature Malleability |
| **Repository** | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/signature-malleability.md) (⭐ 90) |
| **Original Path** | `plugins/scv-scan/skills/scv-scan/references/signature-malleability.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-13 |
| **Updated** | 2026-02-13 |
| **File Hash** | `2fd4c4c60e194ff9...` |

## Description

function claimReward(bytes memory signature, uint256 amount) external {
    // Deduplication by raw signature bytes
    require(!usedSignatures[signature], "already used");

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/signature-malleability.md)*
