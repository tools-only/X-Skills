# Gotchas

| Property | Value |
|----------|-------|
| **Name** | Gotchas |
| **Repository** | [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/pulumi/gotchas.md) (⭐ 392) |
| **Original Path** | `packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/pulumi/gotchas.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-05 |
| **Updated** | 2026-02-05 |
| **File Hash** | `3bc73c7c730f602a...` |

## Description

Problem: Worker fails with "Cannot use import statement outside a module"  
Cause: Pulumi doesn't bundle Worker code  uploads exactly what you provide  
Solution: Build Worker BEFORE Pulumi deploy

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/pulumi/gotchas.md)*
