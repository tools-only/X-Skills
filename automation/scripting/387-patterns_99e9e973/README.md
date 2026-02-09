# Patterns

| Property | Value |
|----------|-------|
| **Name** | Patterns |
| **Repository** | [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/cron-triggers/patterns.md) (🔥 19.8k) |
| **Original Path** | `cli-tool/components/skills/development/cloudflare-deploy/references/cron-triggers/patterns.md` |
| **Category** | automation |
| **Subcategory** | scripting |
| **Tags** | automation |
| **Created** | 2026-02-08 |
| **Updated** | 2026-02-08 |
| **File Hash** | `99e9e973fee96b3a...` |

## Description

typescript
export default {
  async scheduled(controller, env, ctx) {
    const result = await env.DB.prepare(DELETE FROM sessions WHERE expires_at < datetime('now')).run();
    console.log(Deleted ${result.meta.changes} expired sessions);
    ctx.waitUntil(env.DB.prepare("VACUUM").run());
  },
};

**Tags:** `automation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/cron-triggers/patterns.md)*
