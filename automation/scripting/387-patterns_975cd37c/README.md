# Patterns

| Property | Value |
|----------|-------|
| **Name** | Patterns |
| **Repository** | [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/pulumi/patterns.md) (🔥 19.8k) |
| **Original Path** | `cli-tool/components/skills/development/cloudflare-deploy/references/pulumi/patterns.md` |
| **Category** | automation |
| **Subcategory** | scripting |
| **Tags** | automation |
| **Created** | 2026-02-08 |
| **Updated** | 2026-02-08 |
| **File Hash** | `975cd37cd6735d79...` |

## Description

typescript
class WorkerApp extends pulumi.ComponentResource {
    constructor(name: string, args: WorkerAppArgs, opts?) {
        super("custom:cloudflare:WorkerApp", name, {}, opts);
        const defaultOpts = {parent: this};

**Tags:** `automation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/pulumi/patterns.md)*
