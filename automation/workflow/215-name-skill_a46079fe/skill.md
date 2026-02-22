---
name: start-milestone
description: "Begin active work on a GitHub milestone. Args: {milestone-number}. Lists all open issues in the milestone, shows current label state, asks for confirmation, then bulk-transitions status labels from 'status:needs-grooming' to 'status:in-progress' and updates GitHub Projects V2 Status to 'In Progress'. Use after group-items-to-milestone when the team is ready to begin the sprint or release cycle."
argument-hint: '{milestone-number}'
user-invocable: true
---
# Start Milestone

Transition a milestone from planning to active: bulk-update labels and Project board status.

API references: [milestones.md](../gh/references/milestones.md) | [projects-v2.md](../gh/references/projects-v2.md)

## Arguments

`$ARGUMENTS` — milestone number (required).

```text
/start-milestone 3
```

## Workflow

### Step 1: Resolve Milestone

```bash
gh api repos/Jamie-BitFlight/claude_skills/milestones/{number} \
  --jq '[.number, .title, .state, .open_issues, .closed_issues] | @tsv'
```

If milestone is closed, report and stop.

If `open_issues == 0`, warn: "No open issues. Add items first with `/group-items-to-milestone {number}`"

### Step 2: List Issues

```bash
gh issue list -R Jamie-BitFlight/claude_skills \
  --milestone "{title}" \
  --state open \
  --json number,title,labels \
  --jq '.[] | [.number, .title, (.labels | map(.name) | join(","))] | @tsv'
```

### Step 3: Confirm

```text
Start milestone "{title}" — {N} open issues

Issues to transition (status:needs-grooming → status:in-progress):
  #12  SAM: Error Recovery        [priority:p1, status:needs-grooming]
  #8   gitlab-skill: Remove URL   [priority:p1, status:needs-grooming]
  #14  create-backlog-item skill  [priority:p2, status:needs-grooming]

Proceed?
```

Use `AskUserQuestion` with Yes / No options. Stop without changes if No.

### Step 4: Ensure Label Exists

```bash
gh label create "status:in-progress" \
  --color "#0075ca" \
  --description "Work is actively in progress" \
  -R Jamie-BitFlight/claude_skills 2>/dev/null || true
```

### Step 5: Bulk Label Transition

Use the Python automation script (preferred — avoids per-issue bash loops):

```bash
uv run .claude/skills/gh/scripts/github_project_setup.py milestone start \
  --number {number} --repo Jamie-BitFlight/claude_skills
```

Or, to preview changes without applying them:

```bash
uv run .claude/skills/gh/scripts/github_project_setup.py milestone start \
  --number {number} --dry-run
```

The script removes `status:needs-grooming` and adds `status:in-progress` to each
open issue, creates the label if missing, reports per-issue results, and exits
non-zero if any issue fails.

**Fallback** — if the script is unavailable, edit issues individually:

```bash
gh issue edit {number} \
  -R Jamie-BitFlight/claude_skills \
  --add-label "status:in-progress" \
  --remove-label "status:needs-grooming"
```

Log successes and failures; continue on per-issue failures.

### Step 6: Update Project V2

If a GitHub Project exists (`gh project list --owner Jamie-BitFlight`), set Status = `In Progress` for each issue's project item via GraphQL — see [projects-v2.md](../gh/references/projects-v2.md).

### Step 7: Report

```text
Milestone #{N} "{title}" started.

  {count} issues transitioned to status:in-progress
  {failed_count} issues failed (see details above)

Work on individual items:
  /work-backlog-item {title}

Track progress:
  gh issue list -R Jamie-BitFlight/claude_skills --milestone "{title}" --state open
```

## Error Handling

- Milestone not found: list open milestones and stop.
- Milestone already closed: report closed date and stop.
- Label edit fails for specific issue: log and continue; report failures at end.
- User declines confirmation: stop without changes.
- `gh` not installed: run `uv run .claude/skills/gh/scripts/setup_gh.py` first.
