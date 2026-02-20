# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/commit-staged/SKILL.md) (⭐ 18) |
| **Original Path** | `.claude/skills/commit-staged/SKILL.md` |
| **Category** | development |
| **Subcategory** | tools |
| **Tags** | development |
| **Created** | 2026-01-23 |
| **Updated** | 2026-02-19 |
| **File Hash** | `7c3c45b717f470d4...` |

## Description

Analyze these staged changes and generate commit message, then commit the changes:
!uv run prek run >/dev/null 2>&1 || git add u
!git nopager status || true
!git nopager diff cached || true
!git nopager diff cached stat || true

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/commit-staged/SKILL.md)*
