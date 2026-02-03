# 2026 01 12 Missing Cached Tokens

| Property | Value |
|----------|-------|
| **Name** | 2026 01 12 Missing Cached Tokens |
| **Repository** | [alchemiststudiosDOTai/tunacode](https://raw.githubusercontent.com/alchemiststudiosDOTai/tunacode/master/.claude/debug_history/2026-01-12_missing-cached-tokens.md) (⭐ 112) |
| **Original Path** | `.claude/debug_history/2026-01-12_missing-cached-tokens.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-01-12 |
| **Updated** | 2026-01-12 |
| **File Hash** | `58855d814773f316...` |

## Description

The new usage tracker accessed cached_tokens directly on provider usage objects, which caused an AttributeError for providers that omit that field. We now normalize usage objects at the boundary and default missing fields to 0, matching the old behavior.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [alchemiststudiosDOTai/tunacode](https://raw.githubusercontent.com/alchemiststudiosDOTai/tunacode/master/.claude/debug_history/2026-01-12_missing-cached-tokens.md)*
