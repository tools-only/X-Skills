# Api

| Property | Value |
|----------|-------|
| **Name** | Api |
| **Repository** | [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/vectorize/api.md) (⭐ 392) |
| **Original Path** | `packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/vectorize/api.md` |
| **Category** | automation |
| **Subcategory** | scripting |
| **Tags** | automation |
| **Created** | 2026-02-05 |
| **Updated** | 2026-02-05 |
| **File Hash** | `31c2695252bccbf6...` |

## Description

typescript
interface VectorizeVector {
  id: string // Max 64 bytes
  values: number[] // Must match index dimensions
  namespace?: string // Optional partition (max 64 bytes)
  metadata?: Record<string, any> // Max 10 KiB
}

**Tags:** `automation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/vectorize/api.md)*
