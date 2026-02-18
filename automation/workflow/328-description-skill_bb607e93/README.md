# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/work-backlog-item/SKILL.md) (⭐ 18) |
| **Original Path** | `.claude/skills/work-backlog-item/SKILL.md` |
| **Category** | automation |
| **Subcategory** | workflow |
| **Tags** | automation |
| **Created** | 2026-02-14 |
| **Updated** | 2026-02-14 |
| **File Hash** | `bb607e93531b1a53...` |

## Description

Invoke when the user starts work on a backlog item — TRIGGER item title substring provided; reads .claude/BACKLOG.md, auto-grooms if no manifest exists, runs RT-ICA to BLOCK on missing inputs before SAM planning, invokes add-new-feature, then writes Plan reference back to BACKLOG.md. STOPS if item already has a Plan field or RT-ICA returns BLOCKED.

**Tags:** `automation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/work-backlog-item/SKILL.md)*
