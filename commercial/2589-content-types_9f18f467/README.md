# Content Types

| Property | Value |
|----------|-------|
| **Name** | Content Types |
| **Repository** | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/shopify/skills/shopify-content/references/content-types.md) (⭐ 536) |
| **Original Path** | `plugins/shopify/skills/shopify-content/references/content-types.md` |
| **Category** | commercial |
| **Subcategory** | ecommerce |
| **Tags** | commercial |
| **Created** | 2026-02-22 |
| **Updated** | 2026-02-22 |
| **File Hash** | `9f18f467f5977795...` |

## Description

| Operation | Method | Endpoint/Mutation |
||||
| List | GraphQL | { pages(first: 50) { edges { node { id title } } } } |
| Create | GraphQL | pageCreate(page: { ... }) |
| Update | GraphQL | pageUpdate(page: { id, ... }) |
| Delete | GraphQL | pageDelete(id: "...") |

**Tags:** `commercial`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/shopify/skills/shopify-content/references/content-types.md)*
