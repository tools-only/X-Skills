# Headless Setup

| Property | Value |
|----------|-------|
| **Name** | Headless Setup |
| **Repository** | [letta-ai/skills](https://raw.githubusercontent.com/letta-ai/skills/main/tools/obsidian-cli/references/headless-setup.md) (⭐ 49) |
| **Original Path** | `tools/obsidian-cli/references/headless-setup.md` |
| **Category** | daily-assistant |
| **Subcategory** | scheduling |
| **Tags** | daily assistant |
| **Created** | 2026-02-10 |
| **Updated** | 2026-02-10 |
| **File Hash** | `1aa038f73f715acd...` |

## Description

Snap confinement creates separate filesystem namespaces. The CLI process and the running Obsidian instance end up with different SingletonLock/SingletonSocket paths, so CLI can't connect via IPC. The deb package avoids this entirely.

**Tags:** `daily assistant`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [letta-ai/skills](https://raw.githubusercontent.com/letta-ai/skills/main/tools/obsidian-cli/references/headless-setup.md)*
