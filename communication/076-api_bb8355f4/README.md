# Api

| Property | Value |
|----------|-------|
| **Name** | Api |
| **Repository** | [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/email-routing/api.md) (🔥 19.8k) |
| **Original Path** | `cli-tool/components/skills/development/cloudflare-deploy/references/email-routing/api.md` |
| **Category** | communication |
| **Subcategory** | email |
| **Tags** | communication |
| **Created** | 2026-02-08 |
| **Updated** | 2026-02-08 |
| **File Hash** | `bb8355f4ef58ab15...` |

## Description

typescript
interface ExportedHandler<Env = unknown> {
  email?(message: ForwardableEmailMessage, env: Env, ctx: ExecutionContext): void | Promise<void>;
}

**Tags:** `communication`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/email-routing/api.md)*
