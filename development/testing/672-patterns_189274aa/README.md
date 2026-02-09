# Patterns

| Property | Value |
|----------|-------|
| **Name** | Patterns |
| **Repository** | [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/bindings/patterns.md) (🔥 19.8k) |
| **Original Path** | `cli-tool/components/skills/development/cloudflare-deploy/references/bindings/patterns.md` |
| **Category** | development |
| **Subcategory** | testing |
| **Tags** | development |
| **Created** | 2026-02-08 |
| **Updated** | 2026-02-08 |
| **File Hash** | `189274aa2bdf6057...` |

## Description

typescript
// authworker
export default {
  async fetch(request: Request, env: Env) {
    const token = request.headers.get('Authorization');
    return new Response(JSON.stringify({ valid: await validateToken(token) }));
  }
}

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/bindings/patterns.md)*
