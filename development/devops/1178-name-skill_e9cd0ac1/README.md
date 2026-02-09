# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/clerk-pack/skills/clerk-local-dev-loop/SKILL.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/clerk-pack/skills/clerk-local-dev-loop/SKILL.md` |
| **Category** | development |
| **Subcategory** | devops |
| **Tags** | development |
| **Created** | 2026-01-09 |
| **Updated** | 2026-01-09 |
| **File Hash** | `e9cd0ac1a1f1d963...` |

## Description

async function createTestUser() {
  const user = await clerkClient.users.createUser({
    emailAddress: ['test@example.com'],
    password: 'testpassword123',
    firstName: 'Test',
    lastName: 'User'
  })
  console.log('Created test user:', user.id)
}

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/clerk-pack/skills/clerk-local-dev-loop/SKILL.md)*
