# Gotchas

| Property | Value |
|----------|-------|
| **Name** | Gotchas |
| **Repository** | [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/tunnel/gotchas.md) (🔥 19.8k) |
| **Original Path** | `cli-tool/components/skills/development/cloudflare-deploy/references/tunnel/gotchas.md` |
| **Category** | daily-assistant |
| **Subcategory** | notes |
| **Tags** | daily assistant |
| **Created** | 2026-02-08 |
| **Updated** | 2026-02-08 |
| **File Hash** | `267e3fad1f1a328b...` |

## Description

Cause: Tunnel not running or not connected
Solution:
bash
cloudflared tunnel info mytunnel      Check status
ps aux | grep cloudflared              Verify running
journalctl u cloudflared n 100       Check logs

**Tags:** `daily assistant`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/tunnel/gotchas.md)*
