# P1 Gitlab Skill Remove Hardcoded Corporate Url

| Property | Value |
|----------|-------|
| **Name** | P1 Gitlab Skill Remove Hardcoded Corporate Url |
| **Repository** | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/backlog/p1-gitlab-skill-remove-hardcoded-corporate-url.md) (⭐ 20) |
| **Original Path** | `.claude/backlog/p1-gitlab-skill-remove-hardcoded-corporate-url.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-23 |
| **Updated** | 2026-02-27 |
| **File Hash** | `134d0fed9aa7a2c9...` |

## Description

`validate_glfm.py` lines 152-153 hardcode `https://sourcery.assaabloy.net` as the default GitLab instance URL. `gitlab-ci-local-guide.md` line 51 also references this URL. This leaks a corporate internal URL into a public repository. Replace with a generic placeholder (e.g., `https://gitlab.example.com`) or make the URL a required argument with no default.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/backlog/p1-gitlab-skill-remove-hardcoded-corporate-url.md)*
