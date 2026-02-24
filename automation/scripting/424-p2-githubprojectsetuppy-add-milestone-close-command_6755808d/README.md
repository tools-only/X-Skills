# P2 Githubprojectsetuppy Add Milestone Close Command

| Property | Value |
|----------|-------|
| **Name** | P2 Githubprojectsetuppy Add Milestone Close Command |
| **Repository** | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/backlog/p2-githubprojectsetuppy-add-milestone-close-command.md) (⭐ 20) |
| **Original Path** | `.claude/backlog/p2-githubprojectsetuppy-add-milestone-close-command.md` |
| **Category** | automation |
| **Subcategory** | scripting |
| **Tags** | automation |
| **Created** | 2026-02-23 |
| **Updated** | 2026-02-23 |
| **File Hash** | `6755808d811df552...` |

## Description

`github_project_setup.py` now has `milestone start` which bulk-transitions `status:needs-grooming` → `status:in-progress`. Add a symmetric `milestone close` command for the `complete-milestone` skill. The command should: (1) validate the milestone is open, (2) list all still-open issues (warn if any remain), (3) transition open issues from `status:in-progress` → `status:done` or close them, (4) close the milestone itself via `milestone.edit(state=

**Tags:** `automation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/backlog/p2-githubprojectsetuppy-add-milestone-close-command.md)*
