# P2 Fix Pre Existing Linting Errors In Gitlabcontextpy

| Property | Value |
|----------|-------|
| **Name** | P2 Fix Pre Existing Linting Errors In Gitlabcontextpy |
| **Repository** | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/backlog/p2-fix-pre-existing-linting-errors-in-gitlabcontextpy.md) (⭐ 20) |
| **Original Path** | `.claude/backlog/p2-fix-pre-existing-linting-errors-in-gitlabcontextpy.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-24 |
| **Updated** | 2026-02-24 |
| **File Hash** | `5fd8a3330cf5cee8...` |

## Description

plugins/gitlab-skill/skills/gitlab-skill/scripts/gitlab_context.py has 3 linting errors: PLC0415 x2 (deferred imports inside functions) and S607 (partial executable path in subprocess). Fix: move imports to top-level, use shutil.which() for glab path resolution. Same pattern was applied in fetch_gitlab_mr.py migration.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/backlog/p2-fix-pre-existing-linting-errors-in-gitlabcontextpy.md)*
