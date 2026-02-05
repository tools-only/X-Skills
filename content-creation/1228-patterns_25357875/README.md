# Patterns

| Property | Value |
|----------|-------|
| **Name** | Patterns |
| **Repository** | [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/argo-smart-routing/patterns.md) (⭐ 392) |
| **Original Path** | `packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/argo-smart-routing/patterns.md` |
| **Category** | content-creation |
| **Subcategory** | editing |
| **Tags** | content creation |
| **Created** | 2026-02-05 |
| **Updated** | 2026-02-05 |
| **File Hash** | `25357875a98342c2...` |

## Description

typescript
async function enableOptimalPerformance(client: Cloudflare, zoneId: string) {
  await Promise.all([
    client.argo.smartRouting.edit({ zone_id: zoneId, value: 'on' }),
    client.argo.tieredCaching.edit({ zone_id: zoneId, value: 'on' }),
  ])
}

**Tags:** `content creation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/argo-smart-routing/patterns.md)*
