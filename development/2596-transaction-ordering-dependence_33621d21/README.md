# Transaction Ordering Dependence

| Property | Value |
|----------|-------|
| **Name** | Transaction Ordering Dependence |
| **Repository** | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/transaction-ordering-dependence.md) (⭐ 90) |
| **Original Path** | `plugins/scv-scan/skills/scv-scan/references/transaction-ordering-dependence.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-13 |
| **Updated** | 2026-02-13 |
| **File Hash** | `33621d21340d8efd...` |

## Description

// Onchain secret submission
function submitAnswer(bytes32 answer) external {
    // Answer visible in mempool — anyone can copy it
    require(keccak256(abi.encodePacked(answer)) == targetHash);
    _reward(msg.sender);
}

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/transaction-ordering-dependence.md)*
