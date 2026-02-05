# Gotchas

| Property | Value |
|----------|-------|
| **Name** | Gotchas |
| **Repository** | [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/email-workers/gotchas.md) (⭐ 392) |
| **Original Path** | `packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/email-workers/gotchas.md` |
| **Category** | communication |
| **Subcategory** | email |
| **Tags** | communication |
| **Created** | 2026-02-05 |
| **Updated** | 2026-02-05 |
| **File Hash** | `fc57f85947b1f7b9...` |

## Description

typescript
// ❌ WRONG: Stream consumed twice
const email = await PostalMime.parse(await new Response(message.raw).arrayBuffer())
const rawText = await new Response(message.raw).text() // EMPTY!

**Tags:** `communication`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/email-workers/gotchas.md)*
