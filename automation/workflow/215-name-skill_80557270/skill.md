---
name: complete-milestone
description: "Close a completed GitHub milestone. Args: {milestone-number}. Audits open and closed issues, offers to carry forward open items to a new or existing milestone, closes the GitHub milestone, updates Project V2 Status to Done for closed issues, and generates a completion summary. Use when a sprint or release is finished and needs to be officially closed."
argument-hint: '{milestone-number}'
user-invocable: true
---
# Complete Milestone

Verify completion state, handle stragglers, close the milestone, and generate a summary.

API references: [milestones.md](../gh/references/milestones.md) | [projects-v2.md](../gh/references/projects-v2.md)

## Arguments

`$ARGUMENTS` — milestone number (required).

```text
/complete-milestone 3
```

## Workflow

### Step 1: Resolve and Audit

```bash
gh api repos/Jamie-BitFlight/claude_skills/milestones/{number} \
  --jq '[.number, .title, .state, .open_issues, .closed_issues, .due_on] | @tsv'
```

If milestone already closed, report closed date and stop.

Fetch all issues:

```bash
gh issue list -R Jamie-BitFlight/claude_skills \
  --milestone "{title}" --state open \
  --json number,title,labels

gh issue list -R Jamie-BitFlight/claude_skills \
  --milestone "{title}" --state closed \
  --json number,title,labels,closedAt
```

### Step 2: Show State

```text
Milestone #{N}: {title}
Due: {due_on or "not set"}

  Closed: {closed_count} issues ✓
  Open:   {open_count} issues ✗

Closed:
  #12  SAM: Error Recovery      (closed 2026-02-20)
  #8   gitlab-skill: Remove URL (closed 2026-02-19)

Open (incomplete):
  #14  create-backlog-item skill  [status:in-progress]
  #17  start-milestone skill      [status:needs-grooming]
```

If `open_count == 0`: skip to Step 4.

### Step 3: Handle Open Issues

Ask: "What should happen to the {open_count} open issues?"

```text
A) Carry forward to a new milestone (create it now)
B) Carry forward to an existing open milestone
C) Remove from milestone (leave open, unassigned)
D) Close them as incomplete
```

**Option A — new milestone:**

Prompt for new milestone title and optional due date. Create via:

```bash
gh api repos/Jamie-BitFlight/claude_skills/milestones \
  -X POST -f title="{new title}" -f state="open"
```

Reassign open issues:

```bash
gh api repos/Jamie-BitFlight/claude_skills/issues/{N} \
  -X PATCH -F milestone={new_number}
```

**Option B — existing milestone:**

```bash
gh api repos/Jamie-BitFlight/claude_skills/milestones \
  --jq '.[] | select(.state=="open") | [.number, .title] | @tsv'
```

User picks one; reassign open issues to it.

**Option C — unassign:**

```bash
gh api repos/Jamie-BitFlight/claude_skills/issues/{N} \
  -X PATCH -F milestone=null
```

**Option D — close:**

```bash
gh issue close {N} -R Jamie-BitFlight/claude_skills \
  --comment "Closed incomplete as part of milestone #{M} completion."
```

### Step 4: Close Milestone

Use the Python automation script (preferred — handles label transitions and milestone close atomically):

```bash
uv run .claude/skills/gh/scripts/github_project_setup.py milestone close \
  --number {number} --repo Jamie-BitFlight/claude_skills
```

Or, to preview changes without applying them:

```bash
uv run .claude/skills/gh/scripts/github_project_setup.py milestone close \
  --number {number} --dry-run
```

The script transitions all remaining open issues to `status:done`, closes the milestone, and prints a completion summary. It exits non-zero if any issue label update fails.

**Fallback** — if the script is unavailable, close the milestone directly:

```bash
gh api repos/Jamie-BitFlight/claude_skills/milestones/{number} \
  -X PATCH -f state="closed"
```

### Step 5: Update Project V2

For all closed issues, set Project V2 Status = `Done` via GraphQL — see [projects-v2.md](../gh/references/projects-v2.md).

### Step 6: Completion Report

```text
Milestone #{N} "{title}" — CLOSED

  Completed:       {closed_count}/{total_count} ({pct}%)
  Carried forward: {count} → milestone #{new_N} "{new_title}"
  Removed:         {count}
  Closed incomplete: {count}

{if carry_forward_count > 0}
Next steps:
  /group-items-to-milestone {new_N}  (add more items)
  /start-milestone {new_N}           (begin the new sprint)
```

## Error Handling

- Milestone not found: list open milestones and stop.
- Milestone already closed: report and stop.
- Issue reassignment fails: log per-item error, continue with remaining.
- No open milestones for Option B: fall back to Option A automatically.
- `gh` not installed: run `uv run .claude/skills/gh/scripts/setup_gh.py` first.
