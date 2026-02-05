# Patterns

| Property | Value |
|----------|-------|
| **Name** | Patterns |
| **Repository** | [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/pulumi/patterns.md) (⭐ 392) |
| **Original Path** | `packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/pulumi/patterns.md` |
| **Category** | automation |
| **Subcategory** | scripting |
| **Tags** | automation |
| **Created** | 2026-02-05 |
| **Updated** | 2026-02-05 |
| **File Hash** | `56b9a29d6776534f...` |

## Description

typescript
class WorkerApp extends pulumi.ComponentResource {
  constructor(name: string, args: WorkerAppArgs, opts?) {
    super('custom:cloudflare:WorkerApp', name, {}, opts)
    const defaultOpts = { parent: this }

**Tags:** `automation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/pulumi/patterns.md)*
