# Patterns

| Property | Value |
|----------|-------|
| **Name** | Patterns |
| **Repository** | [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/pages/patterns.md) (⭐ 392) |
| **Original Path** | `packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/pages/patterns.md` |
| **Category** | daily-assistant |
| **Subcategory** | tasks |
| **Tags** | daily assistant |
| **Created** | 2026-02-05 |
| **Updated** | 2026-02-05 |
| **File Hash** | `ca26bc18b0a28427...` |

## Description

try {
    const payload = await verifyJWT(authHeader.substring(7), context.env.JWT_SECRET)
    context.data.user = payload
    return context.next()
  } catch (err) {
    return new Response('Invalid token', { status: 401 })
  }
}
export const onRequest = [auth]

**Tags:** `daily assistant`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/pages/patterns.md)*
