# Api

| Property | Value |
|----------|-------|
| **Name** | Api |
| **Repository** | [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/email-routing/api.md) (⭐ 392) |
| **Original Path** | `packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/email-routing/api.md` |
| **Category** | communication |
| **Subcategory** | email |
| **Tags** | communication |
| **Created** | 2026-02-05 |
| **Updated** | 2026-02-05 |
| **File Hash** | `686b9105c873ac87...` |

## Description

typescript
interface ExportedHandler<Env = unknown> {
  email?(message: ForwardableEmailMessage, env: Env, ctx: ExecutionContext): void | Promise<void>
}

**Tags:** `communication`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/email-routing/api.md)*
