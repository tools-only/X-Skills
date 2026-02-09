# Configuration

| Property | Value |
|----------|-------|
| **Name** | Configuration |
| **Repository** | [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/workerd/configuration.md) (🔥 19.8k) |
| **Original Path** | `cli-tool/components/skills/development/cloudflare-deploy/references/workerd/configuration.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-08 |
| **Updated** | 2026-02-08 |
| **File Hash** | `656371e6aa55e315...` |

## Description

const config :Workerd.Config = (
  services = [(name = "main", worker = .mainWorker)],
  sockets = [(name = "http", address = ":8080", http = (), service = "main")]
);

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/workerd/configuration.md)*
