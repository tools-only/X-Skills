# Insufficient Gas Griefing

| Property | Value |
|----------|-------|
| **Name** | Insufficient Gas Griefing |
| **Repository** | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/insufficient-gas-griefing.md) (⭐ 90) |
| **Original Path** | `plugins/scv-scan/skills/scv-scan/references/insufficient-gas-griefing.md` |
| **Category** | daily-assistant |
| **Subcategory** | notes |
| **Tags** | daily assistant |
| **Created** | 2026-02-13 |
| **Updated** | 2026-02-13 |
| **File Hash** | `473c05303de16061...` |

## Description

// Relayer can provide just enough gas for the outer tx to succeed
    // but insufficient gas for the inner call — it silently fails
    (bool success,) = target.call{gas: gasLimit}(data);
    // success is false, but the nonce is already consumed
    // The action is permanently censored
}

**Tags:** `daily assistant`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/insufficient-gas-griefing.md)*
