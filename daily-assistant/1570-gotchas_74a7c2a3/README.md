# Gotchas

| Property | Value |
|----------|-------|
| **Name** | Gotchas |
| **Repository** | [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/durable-objects/gotchas.md) (⭐ 392) |
| **Original Path** | `packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/durable-objects/gotchas.md` |
| **Category** | daily-assistant |
| **Subcategory** | notes |
| **Tags** | daily assistant |
| **Created** | 2026-02-05 |
| **Updated** | 2026-02-05 |
| **File Hash** | `74a7c2a30abfa11b...` |

## Description

Problem: Variables lost after hibernation  
Cause: DO autohibernates when idle; inmemory state not persisted  
Solution: Use ctx.storage for critical data, ws.serializeAttachment() for perconnection metadata

**Tags:** `daily assistant`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/durable-objects/gotchas.md)*
