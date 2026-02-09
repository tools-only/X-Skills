# Gotchas

| Property | Value |
|----------|-------|
| **Name** | Gotchas |
| **Repository** | [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/workers/gotchas.md) (🔥 19.8k) |
| **Original Path** | `cli-tool/components/skills/development/cloudflare-deploy/references/workers/gotchas.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-08 |
| **Updated** | 2026-02-08 |
| **File Hash** | `f9430e9fe0462c57...` |

## Description

Cause: Worker exceeded CPU time limit (10ms standard, 30ms unbound)  
Solution: Use ctx.waitUntil() for background work, offload heavy compute to Durable Objects, or consider Workers AI for ML workloads

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/workers/gotchas.md)*
