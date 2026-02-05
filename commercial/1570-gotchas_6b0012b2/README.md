# Gotchas

| Property | Value |
|----------|-------|
| **Name** | Gotchas |
| **Repository** | [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/d1/gotchas.md) (⭐ 392) |
| **Original Path** | `packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/d1/gotchas.md` |
| **Category** | commercial |
| **Subcategory** | ecommerce |
| **Tags** | commercial |
| **Created** | 2026-02-05 |
| **Updated** | 2026-02-05 |
| **File Hash** | `6b0012b2a24cf478...` |

## Description

Cause: Using string interpolation instead of prepared statements with bind()  
Solution: ALWAYS use prepared statements: env.DB.prepare('SELECT  FROM users WHERE id = ?').bind(userId).all() instead of string interpolation which allows attackers to inject malicious SQL

**Tags:** `commercial`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/d1/gotchas.md)*
