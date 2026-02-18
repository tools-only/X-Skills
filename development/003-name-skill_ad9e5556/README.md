# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CCSBotManager_AddBot_BotNavIgnore/SKILL.md) (⭐ 19) |
| **Original Path** | `.claude/skills/find-CCSBotManager_AddBot_BotNavIgnore/SKILL.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-18 |
| **Updated** | 2026-02-18 |
| **File Hash** | `ad9e5556606aeb19...` |

## Description

Locate the g_pNavMesh nullcheck at the beginning of CCSBotManager_AddBot and generate a patch that converts the conditional jz (jumpifzero) into an unconditional jmp past the earlyreturn, so bots can be added even when no navigation mesh is loaded.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CCSBotManager_AddBot_BotNavIgnore/SKILL.md)*
