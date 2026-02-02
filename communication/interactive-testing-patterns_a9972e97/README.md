# Interactive Testing Patterns

| Property | Value |
|----------|-------|
| **Name** | Interactive Testing Patterns |
| **Repository** | [mrgoonie/claudekit-skills](https://raw.githubusercontent.com/mrgoonie/claudekit-skills/main/.claude/skills/web-testing/references/interactive-testing-patterns.md) (⭐ 1.5k) |
| **Original Path** | `.claude/skills/web-testing/references/interactive-testing-patterns.md` |
| **Category** | communication |
| **Subcategory** | email |
| **Tags** | communication |
| **Created** | 2026-01-19 |
| **Updated** | 2026-01-19 |
| **File Hash** | `a9972e974aa4b848...` |

## Description

javascript
// Text input validation
test('validates email format', async ({ page }) => {
  await page.fill('[name="email"]', 'invalid');
  await page.click('button[type="submit"]');
  await expect(page.locator('.error')).toContainText('Invalid email');
});

**Tags:** `communication`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [mrgoonie/claudekit-skills](https://raw.githubusercontent.com/mrgoonie/claudekit-skills/main/.claude/skills/web-testing/references/interactive-testing-patterns.md)*
