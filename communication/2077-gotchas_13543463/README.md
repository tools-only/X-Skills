# Gotchas

| Property | Value |
|----------|-------|
| **Name** | Gotchas |
| **Repository** | [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/email-workers/gotchas.md) (🔥 19.8k) |
| **Original Path** | `cli-tool/components/skills/development/cloudflare-deploy/references/email-workers/gotchas.md` |
| **Category** | communication |
| **Subcategory** | email |
| **Tags** | communication |
| **Created** | 2026-02-08 |
| **Updated** | 2026-02-08 |
| **File Hash** | `135434634be61d82...` |

## Description

typescript
// ❌ WRONG: Stream consumed twice
const email = await PostalMime.parse(await new Response(message.raw).arrayBuffer());
const rawText = await new Response(message.raw).text(); // EMPTY!

**Tags:** `communication`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/email-workers/gotchas.md)*
