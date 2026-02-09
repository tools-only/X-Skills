# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/linear-pack/skills/linear-performance-tuning/SKILL.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/linear-pack/skills/linear-performance-tuning/SKILL.md` |
| **Category** | automation |
| **Subcategory** | scripting |
| **Tags** | automation |
| **Created** | 2026-01-09 |
| **Updated** | 2026-01-09 |
| **File Hash** | `f9dd4642beae3c4d...` |

## Description

Minimize Field Selection:
typescript
// BAD: Fetching unnecessary fields
const issues = await client.issues();
for (const issue of issues.nodes) {
  // Only using id and title, but fetching everything
  console.log(issue.id, issue.title);
}

**Tags:** `automation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/linear-pack/skills/linear-performance-tuning/SKILL.md)*
