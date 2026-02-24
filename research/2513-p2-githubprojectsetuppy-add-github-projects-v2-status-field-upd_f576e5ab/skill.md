---
name: 'github_project_setup.py: add GitHub Projects V2 status field updates'
description: The `milestone start` and `milestone close` commands transition `status:*` labels on issues but do not update the GitHub Projects V2 `Status` field (the kanban column). These are separate systems — labels are on the issue, the V2 status is a project-level field. Add `project update-status` subcommand to `github_project_setup.py` that uses the GraphQL API to set the `Status` field on project items. This unlocks the full kanban board view in GitHub Projects.
metadata:
  topic: githubprojectsetuppy-add-github-projects-v2-status-field-upd
  source: 'PR #149 follow-up — start-milestone automation (2026-02-22)'
  added: '2026-02-22'
  priority: P2
  type: Feature
  status: open
---
**Suggested location**: `.claude/skills/gh/scripts/github_project_setup.py` — add `project_app` Typer sub-app with `project update-status --project-number N --issue-number N --status "In Progress"` command

**Research first**: Projects V2 GraphQL API — `updateProjectV2ItemFieldValue` mutation. Field ID discovery via `projectV2Fields` query. Reference: `.claude/skills/gh/references/projects-v2.md`.
