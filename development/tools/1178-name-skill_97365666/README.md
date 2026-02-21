# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/development-harness/skills/implementation-manager/SKILL.md) (⭐ 18) |
| **Original Path** | `plugins/development-harness/skills/implementation-manager/SKILL.md` |
| **Category** | development |
| **Subcategory** | tools |
| **Tags** | development |
| **Created** | 2026-02-16 |
| **Updated** | 2026-02-20 |
| **File Hash** | `973656668c7cfd11...` |

## Description

Active task context (if any):
!python3 c "import pathlib, json; context_dir = pathlib.Path('.claude/context'); files = list(context_dir.glob('activetask.json')) if context_dir.exists() else []; print(files[0].read_text() if files else 'No active task')" 2>/dev/null || echo "No active task"

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/development-harness/skills/implementation-manager/SKILL.md)*
