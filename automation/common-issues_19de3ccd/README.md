# Common Issues

| Property | Value |
|----------|-------|
| **Name** | Common Issues |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/code-reviewer/references/common-issues.md) (⭐ 216) |
| **Original Path** | `skills/code-reviewer/references/common-issues.md` |
| **Category** | automation |
| **Subcategory** | scripting |
| **Tags** | automation |
| **Created** | 2025-12-14 |
| **Updated** | 2026-01-29 |
| **File Hash** | `19de3ccd5611e82d...` |

## Description

typescript
// N+1 queries  BAD
const posts = await Post.findAll();
for (const post of posts) {
  post.author = await User.findById(post.authorId); // N queries!
}

**Tags:** `automation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/code-reviewer/references/common-issues.md)*
