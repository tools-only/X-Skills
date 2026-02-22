---
name: group-items-to-milestone
description: "Use when assigning BACKLOG.md items to a GitHub milestone. Args: {milestone-number} [P0|P1|P2|title-filter]. Reads BACKLOG.md, shows items with GitHub Issue status, lets user select which to assign. Creates missing GitHub Issues for selected P0/P1 items, assigns all to the milestone, updates Project V2 Status to Backlog. Use after create-milestone to populate a sprint or release."
argument-hint: '{milestone-number} [P0|P1|P2|title-filter]'
user-invocable: true
---
# Group Items to Milestone

Assign BACKLOG.md items to a GitHub milestone. Bridges BACKLOG.md → GitHub Issues → milestone assignment.

API references: [milestones.md](../gh/references/milestones.md) | [issue-stories.md](../gh/references/issue-stories.md) | [projects-v2.md](../gh/references/projects-v2.md)

## Arguments

- `{milestone-number}` — required
- Optional filter: `P0`, `P1`, `P2`, or title substring to pre-filter the list

```text
/group-items-to-milestone 3
/group-items-to-milestone 3 P1
/group-items-to-milestone 3 github
```

## Workflow

### Step 1: Resolve Milestone

```bash
gh api repos/Jamie-BitFlight/claude_skills/milestones/{number} \
  --jq '[.number, .title, .state] | @tsv'
```

If milestone not found or closed, report and stop.

### Step 2: Load Backlog Items

Read `.claude/BACKLOG.md`. Parse all H3 items from P0, P1, P2, and Ideas sections. Apply any filter.

For each item determine status:

- **Has issue** — `**Issue**: #N` field present → verify state via `gh issue view N --json state`
- **No issue** — P0/P1 item without issue → flagged for creation offer
- **Already in milestone** — issue already assigned to this milestone → shown pre-checked

### Step 3: Present Selection

```text
Milestone #{N}: {title}

P0
  1. [✓] SAM: Error Recovery — Issue #12 (open)
  2. [ ] bash-development: Fix inaccuracies — no issue yet

P1
  3. [✓] gitlab-skill: Remove URL — Issue #8 (open)
  4. [ ] create-backlog-item skill — no issue yet
  5. [~] commitlint verify flag — Issue #5 (already in this milestone)

Legend: [✓] has issue  [ ] needs issue created  [~] already assigned
```

Use `AskUserQuestion`: "Which items to add? (comma-separated numbers, or 'all', or 'P0', 'P1')"

### Step 4: Create Missing Issues

For each selected item with no `**Issue**: #N`:

Build story-format body (Story / Description / Acceptance Criteria / Context). Create issue using the Python script (preferred — handles label creation automatically):

```bash
uv run .claude/skills/gh/scripts/github_project_setup.py issue create \
  --repo Jamie-BitFlight/claude_skills \
  --title "{type}: {title}" \
  --body "{story body}" \
  --priority-label "priority:{p0|p1|p2|idea}" \
  --type-label "type:{feature|bug|refactor|docs|chore}" \
  --milestone {number}
```

Write `**Issue**: #{N}` back to the matching item in `.claude/BACKLOG.md`.

Skip issue creation for P2/Ideas items — assign by milestone number only if they already have an issue.

### Step 5: Assign Existing Issues

For selected items that already have issues but are not yet in this milestone:

```bash
gh api repos/Jamie-BitFlight/claude_skills/issues/{issue_number} \
  -X PATCH -F milestone={milestone_number}
```

### Step 6: Update Project V2 Status

If a GitHub Project exists (`gh project list --owner Jamie-BitFlight`), set Status = `Backlog` for each newly assigned item via GraphQL — see [projects-v2.md](../gh/references/projects-v2.md).

### Step 7: Report

```text
Milestone #{N}: {title}

Assigned {count} items:
  Issue #12: SAM: Error Recovery (existing issue)
  Issue #14: create-backlog-item skill (new issue created)
  Issue #5:  commitlint verify flag (already assigned — skipped)

BACKLOG.md updated with {created_count} new issue numbers.

Next step: /start-milestone {number}
```

## Error Handling

- Milestone not found: list open milestones and stop.
- Issue creation fails: log error per item, continue with remaining.
- `gh` not installed: run `uv run .claude/skills/gh/scripts/setup_gh.py` first.
- No items match filter: report and show available sections.
- Label not found: create it on the fly with `gh label create` then retry.
