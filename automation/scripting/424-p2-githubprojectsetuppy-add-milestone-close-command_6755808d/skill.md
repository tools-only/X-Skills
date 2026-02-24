---
name: 'github_project_setup.py: add `milestone close` command'
description: "`github_project_setup.py` now has `milestone start` which bulk-transitions `status:needs-grooming` → `status:in-progress`. Add a symmetric `milestone close` command for the `complete-milestone` skill. The command should: (1) validate the milestone is open, (2) list all still-open issues (warn if any remain), (3) transition open issues from `status:in-progress` → `status:done` or close them, (4) close the milestone itself via `milestone.edit(state='closed')`, (5) print a completion summary. Updat"
metadata:
  topic: githubprojectsetuppy-add-milestone-close-command
  source: 'PR #149 follow-up — start-milestone automation (2026-02-22)'
  added: '2026-02-22'
  priority: P2
  type: Feature
  status: done
---
**Suggested location**: `.claude/skills/gh/scripts/github_project_setup.py` — add `milestone close` subcommand under `milestone_app`; update `.claude/skills/complete-milestone/SKILL.md`
