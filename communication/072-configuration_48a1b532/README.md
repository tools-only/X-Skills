# Configuration

| Property | Value |
|----------|-------|
| **Name** | Configuration |
| **Repository** | [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/email-workers/configuration.md) (⭐ 392) |
| **Original Path** | `packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/email-workers/configuration.md` |
| **Category** | communication |
| **Subcategory** | email |
| **Tags** | communication |
| **Created** | 2026-02-05 |
| **Updated** | 2026-02-05 |
| **File Hash** | `48a1b532c9d120d6...` |

## Description

typescript
interface Env {
  EMAIL: SendEmail
  ARCHIVE: KVNamespace
  ATTACHMENTS: R2Bucket
  WEBHOOK_URL: string
}

**Tags:** `communication`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/email-workers/configuration.md)*
