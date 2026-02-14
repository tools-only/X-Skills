# Unsecure Signatures

| Property | Value |
|----------|-------|
| **Name** | Unsecure Signatures |
| **Repository** | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/unsecure-signatures.md) (⭐ 90) |
| **Original Path** | `plugins/scv-scan/skills/scv-scan/references/unsecure-signatures.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-13 |
| **Updated** | 2026-02-13 |
| **File Hash** | `44d435debb3021ac...` |

## Description

// 4. Raw ecrecover — no null address check
    address recovered = ecrecover(hash, v, r, s);
    // 5. No svalue malleability check
    require(recovered == signer);

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/unsecure-signatures.md)*
