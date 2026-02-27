# GitHub Integration — Detailed Steps and Examples

## Step 2.5: GitHub Issue Sync

<github_sync>

After extracting item fields (Step 2), check for an existing linked issue:

1. Search the matched item for `**Issue**: #N` field.
2. If found: issue already linked. Run:

   ```bash
   gh issue view N -R Jamie-BitFlight/claude_skills --json number,title,state,labels
   ```

   Report the issue state. If open, proceed. If closed, warn the user before re-opening planning.

3. If not found AND priority is P0 or P1: offer to create a GitHub Issue:

   ```text
   This P1 item has no linked GitHub issue. Create one? (yes/no)
   ```

   If yes, proceed to 2.5a.
   If no, skip GitHub sync; the per-item file remains the only local record.

4. If not found AND priority is P2 or Ideas: do not prompt; skip GitHub sync silently.

</github_sync>

## Step 2.5a: Create GitHub Issue

Invoke the backlog script:

```bash
uv run .claude/skills/backlog/scripts/backlog.py update "{title}" --create-issue -R Jamie-BitFlight/claude_skills
```

The script creates the issue and writes `issue: '#N'` back to the per-item file frontmatter.

## Step 2.7: Set In-Progress Label

Invoke the backlog script:

```bash
uv run .claude/skills/backlog/scripts/backlog.py update "{title}" --status in-progress -R Jamie-BitFlight/claude_skills
```

If the item is in a milestone with other issues, also run `milestone start`:

```bash
uv run .claude/skills/gh/scripts/github_project_setup.py milestone start \
  --number {milestone_number} --repo Jamie-BitFlight/claude_skills
```

## Step 9: Close — backlog script

Invoke the backlog script (updates per-item file + gh issue close):

```bash
uv run .claude/skills/backlog/scripts/backlog.py close "{title}" --plan "{plan path}" --checklist-pass -R Jamie-BitFlight/claude_skills
```

The script updates the per-item file status and closes the GitHub issue.

## setup-github Command

**Trigger:** `$ARGUMENTS` is exactly `setup-github`.

<setup_github>

1. Run label taxonomy setup:

   ```bash
   uv run .claude/skills/gh/scripts/github_project_setup.py labels \
     --repo Jamie-BitFlight/claude_skills
   ```

2. Check for existing milestones:

   ```bash
   gh api repos/Jamie-BitFlight/claude_skills/milestones
   ```

   If none exist, create the first milestone:

   ```bash
   gh api repos/Jamie-BitFlight/claude_skills/milestones \
     -X POST \
     -f title="v1.0 — Skills Foundation" \
     -f description="Initial stable milestone for claude_skills skills and plugins" \
     -f due_on="2026-03-31T00:00:00Z"
   ```

3. Check for existing projects:

   ```bash
   gh project list --owner Jamie-BitFlight
   ```

   If none exist, prompt: "Create GitHub Project 'claude_skills Backlog'? (yes/no)"
   If yes:

   ```bash
   gh project create --owner Jamie-BitFlight --title "claude_skills Backlog"
   gh project link 1 --owner Jamie-BitFlight --repo Jamie-BitFlight/claude_skills
   ```

4. Report setup summary:

   ```text
   GitHub setup complete:
   - Labels: N created
   - Milestone: #1 "v1.0 — Skills Foundation"
   - Project: #1 "claude_skills Backlog" (linked to repo)

   Next steps:
   - Add custom fields: .claude/skills/gh/references/projects-v2.md
   - Import existing backlog: /work-backlog-item <title> for each P0/P1 item
   ```

</setup_github>

---

## Example Sessions

### GitHub Issue Creation Flow

```text
> /work-backlog-item clang-format YAML frontmatter

Found: "clang-format: Fix broken YAML frontmatter" (P1)
No grooming manifest. Running groom-backlog-item first...

RT-ICA: APPROVED
This P1 item has no linked GitHub issue. Create one? yes

Creating issue...
  Title: fix: correct YAML frontmatter in clang-format SKILL.md
  Labels: priority:p1, type:bug, status:needs-grooming

Created: #59 — https://github.com/Jamie-BitFlight/claude_skills/issues/59
Updated per-item file: issue: '#59'

Setting status:in-progress...
  ✓ Added status:in-progress, removed status:needs-grooming

Composing feature request...
Invoking /python3-development:add-new-feature...

[SAM phases run]

Updated per-item file with Plan: plan/tasks-3-clang-format-yaml.md
GitHub issue #59 — Plan added to issue body.

Next steps:
- To execute:      /python3-development:implement-feature clang-format-yaml
- To close when done: /work-backlog-item close clang-format
```

### GitHub Setup Session

```text
> /work-backlog-item setup-github

Setting up GitHub project infrastructure...

1. Creating label taxonomy...
   created: priority:p0, priority:p1, priority:p2, priority:idea
   created: type:feature, type:bug, type:refactor, type:docs, type:chore
   created: status:in-progress, status:blocked, status:needs-grooming, status:needs-review
   Labels: 13 created, 0 skipped

2. Creating milestone: v1.0 — Skills Foundation
   Created milestone #1

3. Create GitHub Project 'claude_skills Backlog'? yes
   Created project #1, linked to repo

GitHub setup complete.
```

---

## Field Mapping Reference

```text
.claude/backlog/    →  GitHub Issue
  metadata.priority →  priority:* label
  Description       →  Issue body (story format)
  metadata.status   →  status:* label
  metadata.plan     →  Issue body Notes
  metadata.issue    ←  written back after creation
  metadata.status   →  Issue closed
```

See [issue-stories.md](.claude/skills/gh/references/issue-stories.md) for the full body template and lifecycle.
