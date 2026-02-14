# Reentrancy

| Property | Value |
|----------|-------|
| **Name** | Reentrancy |
| **Repository** | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/reentrancy.md) (⭐ 90) |
| **Original Path** | `plugins/scv-scan/skills/scv-scan/references/reentrancy.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-13 |
| **Updated** | 2026-02-13 |
| **File Hash** | `44db9f255801f719...` |

## Description

// Hidden external calls that trigger callbacks:
// ERC721._safeMint() > onERC721Received()
// ERC1155.safeTransferFrom() > onERC1155Received()
// ERC777 token transfers > tokensReceived() hook

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/reentrancy.md)*
