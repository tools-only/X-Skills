# Gotchas

| Property | Value |
|----------|-------|
| **Name** | Gotchas |
| **Repository** | [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/workers/gotchas.md) (⭐ 392) |
| **Original Path** | `packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/workers/gotchas.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-05 |
| **Updated** | 2026-02-05 |
| **File Hash** | `c82806cd323584c2...` |

## Description

Cause: Worker exceeded CPU time limit (10ms standard, 30ms unbound)  
Solution: Use ctx.waitUntil() for background work, offload heavy compute to Durable Objects, or consider Workers AI for ML workloads

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/workers/gotchas.md)*
