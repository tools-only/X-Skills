# Gotchas

| Property | Value |
|----------|-------|
| **Name** | Gotchas |
| **Repository** | [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/workerd/gotchas.md) (🔥 19.8k) |
| **Original Path** | `cli-tool/components/skills/development/cloudflare-deploy/references/workerd/gotchas.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-08 |
| **Updated** | 2026-02-08 |
| **File Hash** | `ed3107537e54d440...` |

## Description

✅ Correct:
capnp
const worker :Workerd.Worker = (
  serviceWorkerScript = embed "worker.js",
  compatibilityDate = "20240115"   Always set!
)

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/workerd/gotchas.md)*
