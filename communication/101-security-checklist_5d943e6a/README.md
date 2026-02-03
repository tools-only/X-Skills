# z.string().min(1).max(100),

| Property | Value |
|----------|-------|
| **Name** | z.string().min(1).max(100), |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/fullstack-guardian/references/security-checklist.md) (⭐ 216) |
| **Original Path** | `skills/fullstack-guardian/references/security-checklist.md` |
| **Category** | communication |
| **Subcategory** | email |
| **Tags** | communication |
| **Created** | 2025-12-14 |
| **Updated** | 2026-01-29 |
| **File Hash** | `5d943e6a3a8f3c34...` |

## Description

typescript
// NestJS Guard
@UseGuards(JwtAuthGuard)
@Get('profile')
async getProfile(@CurrentUser() user: User) {
  return this.userService.findById(user.id);
}

**Tags:** `communication`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/fullstack-guardian/references/security-checklist.md)*
