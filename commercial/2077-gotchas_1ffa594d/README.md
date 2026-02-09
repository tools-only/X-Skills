# Gotchas

| Property | Value |
|----------|-------|
| **Name** | Gotchas |
| **Repository** | [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/d1/gotchas.md) (🔥 19.8k) |
| **Original Path** | `cli-tool/components/skills/development/cloudflare-deploy/references/d1/gotchas.md` |
| **Category** | commercial |
| **Subcategory** | ecommerce |
| **Tags** | commercial |
| **Created** | 2026-02-08 |
| **Updated** | 2026-02-08 |
| **File Hash** | `1ffa594d0bdc7e02...` |

## Description

Cause: Using string interpolation instead of prepared statements with bind()  
Solution: ALWAYS use prepared statements: env.DB.prepare('SELECT  FROM users WHERE id = ?').bind(userId).all() instead of string interpolation which allows attackers to inject malicious SQL

**Tags:** `commercial`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/d1/gotchas.md)*
