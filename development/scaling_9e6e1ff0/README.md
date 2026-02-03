# socket.io-stream

| Property | Value |
|----------|-------|
| **Name** | socket.io-stream |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/websocket-engineer/references/scaling.md) (⭐ 216) |
| **Original Path** | `skills/websocket-engineer/references/scaling.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2025-12-15 |
| **Updated** | 2025-12-15 |
| **File Hash** | `9e6e1ff026caceac...` |

## Description

┌─────────────┐
│Load Balancer│ (nginx/HAProxy with sticky sessions)
└──────┬──────┘
       │
   ┌───┴───┐
   │       │
┌──▼──┐ ┌──▼──┐
│WS 1│ │WS 2│ ... (Socket.IO servers)
└──┬──┘ └──┬──┘
   │       │
   └───┬───┘
       │
   ┌───▼───┐
   │ Redis │ (Pub/Sub adapter)
   └───────┘

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/websocket-engineer/references/scaling.md)*
