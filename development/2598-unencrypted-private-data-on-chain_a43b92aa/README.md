# Unencrypted Private Data On Chain

| Property | Value |
|----------|-------|
| **Name** | Unencrypted Private Data On Chain |
| **Repository** | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/unencrypted-private-data-on-chain.md) (⭐ 90) |
| **Original Path** | `plugins/scv-scan/skills/scv-scan/references/unencrypted-private-data-on-chain.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-13 |
| **Updated** | 2026-02-13 |
| **File Hash** | `a43b92aafdbe28d9...` |

## Description

constructor(bytes32 _answer, string memory _pwd) {
        secretAnswer = _answer; // Visible in deployment tx calldata
        password = _pwd;        // Readable from storage slot
    }

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/unencrypted-private-data-on-chain.md)*
