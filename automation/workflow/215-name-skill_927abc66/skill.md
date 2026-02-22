---
name: create-milestone
description: "Create a GitHub milestone for Jamie-BitFlight/claude_skills. No args: guided intake (title, due date, description). With 'quick {title}': minimal prompts. Creates milestone via gh api REST, checks for duplicates. Returns milestone number for use with group-items-to-milestone. Use when starting a new sprint, release, or theme grouping of backlog items."
argument-hint: '[quick {title}]'
user-invocable: true
---
# Create Milestone

Create a GitHub milestone on `Jamie-BitFlight/claude_skills` and return its number for downstream use.

REST API reference: [milestones.md](../gh/references/milestones.md)

## Arguments

- **Empty** — guided intake via `AskUserQuestion`
- **`quick {title}`** — use remainder as title, ask only for due date and description

## Workflow

### Step 1: Collect Fields

**Guided mode** (no args) — ask in sequence:

```text
Q1: Milestone title?  (e.g. "v1.1 — Milestone Workflow", "2026-Q1 Grooming")
Q2: Due date? (YYYY-MM-DD, or skip)
Q3: Description? (one sentence, or skip)
```

**Quick mode** (`quick {title}`) — ask only Q2 and Q3.

Title is required. Due date and description are optional.

### Step 2: Duplicate Check

```bash
gh api repos/Jamie-BitFlight/claude_skills/milestones \
  --jq '.[] | select(.state=="open") | .title'
```

If an open milestone with the same title already exists, report it and ask: "Use existing or create new?" via `AskUserQuestion`.

If user chooses existing: print its number and stop.

### Step 3: Create Milestone

Use the Python script (preferred — returns structured output):

```bash
uv run .claude/skills/gh/scripts/github_project_setup.py milestone create \
  --repo Jamie-BitFlight/claude_skills \
  --title "{title}" \
  --description "{description}" \
  --due "{YYYY-MM-DD}"
```

Omit `--due` if not provided. Omit `--description` if not provided.

Capture the milestone number from the output line `Created milestone #{number}: …`.

### Step 4: Confirm

```text
Milestone created.

  Title:   {title}
  Number:  #{number}
  Due:     {due date or "not set"}
  URL:     https://github.com/Jamie-BitFlight/claude_skills/milestone/{number}

Next steps:
  Assign items:  /group-items-to-milestone {number}
  Start work:    /start-milestone {number}
```

## Error Handling

- `gh` not installed: run `uv run .claude/skills/gh/scripts/setup_gh.py` first.
- `GITHUB_TOKEN` missing: report and stop.
- Duplicate found and user picks existing: print existing milestone number and next-step commands, stop.
- API error: print full response and stop.
