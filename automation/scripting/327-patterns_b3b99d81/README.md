# Patterns

| Property | Value |
|----------|-------|
| **Name** | Patterns |
| **Repository** | [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/do-storage/patterns.md) (⭐ 392) |
| **Original Path** | `packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/do-storage/patterns.md` |
| **Category** | automation |
| **Subcategory** | scripting |
| **Tags** | automation |
| **Created** | 2026-02-05 |
| **Updated** | 2026-02-05 |
| **File Hash** | `b3b99d810eb3507b...` |

## Description

typescript
export class MyDurableObject extends DurableObject {
  constructor(ctx: DurableObjectState, env: Env) {
    super(ctx, env)
    this.sql = ctx.storage.sql

**Tags:** `automation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/do-storage/patterns.md)*
