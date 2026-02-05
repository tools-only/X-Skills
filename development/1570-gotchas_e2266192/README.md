# Gotchas

| Property | Value |
|----------|-------|
| **Name** | Gotchas |
| **Repository** | [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/pages/gotchas.md) (⭐ 392) |
| **Original Path** | `packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/pages/gotchas.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-05 |
| **Updated** | 2026-02-05 |
| **File Hash** | `e2266192c4a1d596...` |

## Description

Problem: Function endpoints return 404 or don't execute  
Causes: _routes.json excludes path; wrong file extension (.jsx/.tsx); Functions dir not at output root  
Solution: Check _routes.json, rename to .ts/.js, verify build output structure

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/pages/gotchas.md)*
