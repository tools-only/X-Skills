# Implementation

| Property | Value |
|----------|-------|
| **Name** | Implementation |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/automation/mattyp-changelog/skills/changelog-orchestrator/references/implementation.md) (⭐ 1.3k) |
| **Original Path** | `plugins/automation/mattyp-changelog/skills/changelog-orchestrator/references/implementation.md` |
| **Category** | content-creation |
| **Subcategory** | editing |
| **Tags** | content creation |
| **Created** | 2026-01-06 |
| **Updated** | 2026-01-06 |
| **File Hash** | `b86076af8e16a441...` |

## Description

1. Read .changelogconfig.json from the repo root.
2. Validate it with {baseDir}/scripts/validate_config.py.
3. Decide date range:
    Weekly mode: today minus 7 days → today
    Custom mode: use provided start_date/end_date

**Tags:** `content creation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/automation/mattyp-changelog/skills/changelog-orchestrator/references/implementation.md)*
