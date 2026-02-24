---
name: create-backlog-item
description: "Create a new item in .claude/BACKLOG.md. No args: guided intake via AskUserQuestion. With 'quick {title}': inline capture from user description. With '--auto {title}': fully autonomous — derives all fields from research files and task context, skips all AskUserQuestion calls, logs decisions made. Writes item to correct priority section, updates frontmatter counts, offers GitHub Issue creation for P0/P1. Enforces required fields (title, priority, description) and duplicate detection before writing."
argument-hint: '[quick {title} | --auto {title} | <empty for guided intake>]'
user-invocable: true
---
# Create Backlog Item

Capture a new backlog item and append it to `.claude/BACKLOG.md` in the correct priority section.

## Arguments

`$0` selects the operating mode:

| `$0` value | Mode | Remaining args |
|---|---|---|
| (empty) | Guided intake via `AskUserQuestion` | — |
| `quick` | Fast entry — skip title question | `$1`+ = title |
| `--auto` | Fully autonomous — no interactive prompts | `$1`+ = title |

```text
/create-backlog-item                                   # guided intake
/create-backlog-item quick GitHub milestone skill      # quick entry
/create-backlog-item --auto vercel skills npm package  # autonomous (agent use)
```

## Workflow

### Step 1: Collect Item Fields

**If `$0` is `--auto`:**

Title = `$1` onward (all remaining words joined). Do not call `AskUserQuestion`. Instead:

1. Search `research/` recursively for any file whose name or content matches the title (case-insensitive). Read the best match.
2. Search `.claude/BACKLOG.md` for related items to understand existing priority patterns.
3. Derive all fields from the research file, task description, and available context:
   - **Title**: from `$1` onward
   - **Priority**: infer from description urgency keywords (`critical`, `required`, `must` → P1; `nice to have`, `optional` → P2; default P1)
   - **Description**: summarize from research file overview + problem statement
   - **Source**: `"Agent task — auto-derived from research/{filename}"`
   - **Type**: infer from description (`install`, `integrate`, `add` → Feature; default Feature)
4. Log every decision:

```text
[AUTO] Title: {title} — from $1 onward
[AUTO] Priority: P1 — inferred from description (no urgency keywords found, defaulting P1)
[AUTO] Description: derived from research/skill-generation-tools/vercel-labs-skills.md
[AUTO] Source: Agent task — auto-derived from research/skill-generation-tools/vercel-labs-skills.md
[AUTO] Type: Feature — inferred from "integrate" keyword
```

Proceed to Step 2 (validate). Skip Step 7 (GitHub issue) — auto mode does not create GitHub issues unless `--create-issue` is also passed.

**If `$0` is empty (guided intake):**

Use `AskUserQuestion` with two questions:

```text
Question 1: "What is the title of this backlog item?"
  header: "Title"
  options: (free text — no options; user will type)

Question 2: "What priority?"
  header: "Priority"
  options:
    - label: "P0 — Must have"
    - label: "P1 — Should have"
    - label: "P2 — Nice to have"
    - label: "Idea — Not yet sized"
```

Then ask:

```text
Question 3: "Describe the item. Include: what the problem is, what success looks like, and any known constraints or research questions."
  header: "Description"
```

Then ask:

```text
Question 4: "What is the source or trigger for this item? (e.g., code review, user request, session observation — or press Enter to skip)"
  header: "Source"
  options:
    - label: "Code review"
    - label: "User request"
    - label: "Session observation"
    - label: "Research finding"
    - label: "Skip"
```

Then ask:

```text
Question 5: "What type of work is this?"
  header: "Type"
  options:
    - label: "Feature"
    - label: "Bug"
    - label: "Refactor"
    - label: "Docs"
    - label: "Chore"
```

**If `$0` is `quick`:**

Title = `$1` onward. Ask only:

- Priority (Question 2 above)
- Description (Question 3 above)

Source defaults to "Session observation". Type defaults to "Feature".

### Step 2: Validate Inputs

Required fields: `title`, `priority`, `description`.

- If any required field is missing or empty: report which field is missing and stop.
- `title` must be non-empty after trimming whitespace.
- `description` must be non-empty after trimming whitespace.

### Step 3: Duplicate Detection

Read `.claude/BACKLOG.md`. Search all H3 headings for case-insensitive overlap with `title`.

If a match is found within edit distance ≤ 2 tokens (same first 3 words), report:

```text
Possible duplicate: "{existing title}" already exists in {section}.
Proceed anyway? (y/n)
```

Use `AskUserQuestion` with Yes / No options. If No: stop.

**In `--auto` mode**: if a duplicate is found, log `[AUTO] STOP — duplicate detected: "{existing title}" in {section}` and stop without writing. Do not ask.

### Step 4: Compose Item Block

Format today's date as `YYYY-MM-DD` (use system date).

```text
### {title}

**Source**: {source, or "Not specified" if skipped}
**Added**: {YYYY-MM-DD}
**Priority**: {P0|P1|P2|Idea}
**Type**: {type}
**Description**: {description}
```

If research questions were embedded in the description (lines starting with `?` or `Research:`), extract them into a separate `**Research first**:` field:

```text
**Research first**: {extracted questions}
```

### Step 5: Write to BACKLOG.md via backlog script

Build the command. Base:

```bash
uv run .claude/skills/backlog/scripts/backlog.py add \
  --title "{title}" \
  --priority "{priority}" \
  --description "{description}" \
  --source "{source}" \
  --type "{type}" \
  -R Jamie-BitFlight/claude_skills
```

- If research_first is non-empty: append `--research-first "{research_first}"`

**GitHub Issue creation:**

- If priority is P0 or P1 and (guided/quick mode with user said Yes, or `--auto` with `--create-issue` passed): do NOT add `--no-create-issue` (script creates issue by default).
- If priority is P2 or Idea: add `--no-create-issue`.
- If priority is P0 or P1 and user said No (skip): add `--no-create-issue`.
- If `$0` is `--auto` and user did not pass `--create-issue`: add `--no-create-issue`.

### Step 6: Confirm Write

The script outputs the confirmation. Report to user:

```text
Backlog item created.

  Title:    {title}
  Priority: {priority}
  Section:  {section heading}
  Added:    {date}

Next steps:
  Groom:  /groom-backlog-item {title}
  Work:   /work-backlog-item {title}
```

## Error Handling

- Missing required field: report field name, stop.
- Duplicate detected and user says No: stop without writing.
- backlog script fails: report error, stop.
- GITHUB_TOKEN not set (for P0/P1 issue creation): script reports; item still written to BACKLOG.md.

## Completion Criteria

- backlog add invoked successfully
- Item written to correct BACKLOG.md section (script handles)
- Frontmatter counts updated (script handles)
- GitHub Issue created and `**Issue**: #N` written back (P0/P1 only, if --create-issue; script handles)
- Next-step commands shown to user
