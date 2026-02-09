# Gotchas

| Property | Value |
|----------|-------|
| **Name** | Gotchas |
| **Repository** | [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/durable-objects/gotchas.md) (🔥 19.8k) |
| **Original Path** | `cli-tool/components/skills/development/cloudflare-deploy/references/durable-objects/gotchas.md` |
| **Category** | daily-assistant |
| **Subcategory** | notes |
| **Tags** | daily assistant |
| **Created** | 2026-02-08 |
| **Updated** | 2026-02-08 |
| **File Hash** | `b7156d2d5415f367...` |

## Description

Problem: Variables lost after hibernation  
Cause: DO autohibernates when idle; inmemory state not persisted  
Solution: Use ctx.storage for critical data, ws.serializeAttachment() for perconnection metadata

**Tags:** `daily assistant`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/durable-objects/gotchas.md)*
