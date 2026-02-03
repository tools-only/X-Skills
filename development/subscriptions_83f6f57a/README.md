# GraphQL Subscriptions

| Property | Value |
|----------|-------|
| **Name** | GraphQL Subscriptions |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/graphql-architect/references/subscriptions.md) (⭐ 216) |
| **Original Path** | `skills/graphql-architect/references/subscriptions.md` |
| **Category** | data-analysis |
| **Subcategory** | processing |
| **Tags** | data analysis |
| **Created** | 2025-12-15 |
| **Updated** | 2026-01-29 |
| **File Hash** | `83f6f57a6b8692da...` |

## Description

typescript
// schema.graphql
type Subscription {
  postCreated: Post!
  postUpdated(id: ID!): Post!
  commentAdded(postId: ID!): Comment!
  userOnline: User!
}

**Tags:** `data analysis`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/graphql-architect/references/subscriptions.md)*
