# Patterns

| Property | Value |
|----------|-------|
| **Name** | Patterns |
| **Repository** | [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/cron-triggers/patterns.md) (⭐ 392) |
| **Original Path** | `packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/cron-triggers/patterns.md` |
| **Category** | automation |
| **Subcategory** | scripting |
| **Tags** | automation |
| **Created** | 2026-02-05 |
| **Updated** | 2026-02-05 |
| **File Hash** | `b807ff7af5a4a7dd...` |

## Description

typescript
export default {
  async scheduled(controller, env, ctx) {
    const result = await env.DB.prepare(DELETE FROM sessions WHERE expires_at < datetime('now')).run()
    console.log(Deleted ${result.meta.changes} expired sessions)
    ctx.waitUntil(env.DB.prepare('VACUUM').run())
  },
}

**Tags:** `automation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/cron-triggers/patterns.md)*
