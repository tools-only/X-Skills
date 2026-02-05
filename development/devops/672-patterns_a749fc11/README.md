# Patterns

| Property | Value |
|----------|-------|
| **Name** | Patterns |
| **Repository** | [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/tunnel/patterns.md) (⭐ 392) |
| **Original Path** | `packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/tunnel/patterns.md` |
| **Category** | development |
| **Subcategory** | devops |
| **Tags** | development |
| **Created** | 2026-02-05 |
| **Updated** | 2026-02-05 |
| **File Hash** | `a749fc115b99b3d8...` |

## Description

yaml
services:
  cloudflared:
    image: cloudflare/cloudflared:latest
    command: tunnel noautoupdate run token ${TUNNEL_TOKEN}
    restart: unlessstopped

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/tunnel/patterns.md)*
