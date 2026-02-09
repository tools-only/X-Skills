# Gotchas

| Property | Value |
|----------|-------|
| **Name** | Gotchas |
| **Repository** | [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/pages-functions/gotchas.md) (🔥 19.8k) |
| **Original Path** | `cli-tool/components/skills/development/cloudflare-deploy/references/pages-functions/gotchas.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-08 |
| **Updated** | 2026-02-08 |
| **File Hash** | `c623392e07ba0b90...` |

## Description

Problem: ctx.env.MY_BINDING shows type error  
Cause: No type definition for Env  
Solution: Run npx wrangler types or manually define:
typescript
interface Env { MY_BINDING: KVNamespace; }
export const onRequest: PagesFunction<Env> = async (ctx) => { / ... / };

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/pages-functions/gotchas.md)*
