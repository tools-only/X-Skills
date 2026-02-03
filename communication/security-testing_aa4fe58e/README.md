# Security Testing

| Property | Value |
|----------|-------|
| **Name** | Security Testing |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/test-master/references/security-testing.md) (⭐ 216) |
| **Original Path** | `skills/test-master/references/security-testing.md` |
| **Category** | communication |
| **Subcategory** | email |
| **Tags** | communication |
| **Created** | 2025-12-14 |
| **Updated** | 2026-01-29 |
| **File Hash** | `aa4fe58e2e1341fc...` |

## Description

typescript
describe('Authentication Security', () => {
  it('rejects invalid credentials', async () => {
    await request(app)
      .post('/api/login')
      .send({ email: 'user@test.com', password: 'wrong' })
      .expect(401);
  });

**Tags:** `communication`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/test-master/references/security-testing.md)*
