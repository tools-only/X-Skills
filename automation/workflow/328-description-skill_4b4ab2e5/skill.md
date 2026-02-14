---
description: Use when starting work on a backlog item. Bridges BACKLOG.md to SAM planning pipeline — reads item, auto-grooms if needed, composes feature request, invokes add-new-feature, then updates backlog with plan reference.
argument-hint: <item-title-substring>
user-invocable: true
---

# Work Backlog Item

Take a backlog item from `.claude/BACKLOG.md` and feed it into the SAM planning pipeline (`/python3-development:add-new-feature`).

## Arguments

`$ARGUMENTS` is a case-insensitive title substring matching a backlog item.

Examples:

```text
/work-backlog-item Error Recovery
/work-backlog-item regex false positive
/work-backlog-item validate-task-file
```

## Workflow

### Step 1: Find the Backlog Item

Read `.claude/BACKLOG.md`. Search for an H3 heading (`### ...`) whose text contains `$ARGUMENTS` (case-insensitive).

If zero matches: report "No backlog item found matching: {$ARGUMENTS}" and stop.
If multiple matches: list all matches and ask the user to pick one.

Record which priority section (P0, P1, P2, Ideas) the item belongs to.

### Step 2: Extract Item Fields

From the matched item, extract these fields (all are **bold-key**: value format on separate lines):

| Field | Key in BACKLOG.md | Required |
|-------|-------------------|----------|
| title | The H3 heading text | Yes |
| source | `**Source**:` | No |
| added | `**Added**:` | No |
| description | `**Description**:` | Yes |
| research_first | `**Research first**:` | No |
| suggested_location | `**Suggested location**:` | No |
| plan | `**Plan**:` | No |

If the item already has a `**Plan**:` field, report: "This item already has a plan at {path}. Use `/python3-development:implement-feature {path}` to execute it." and stop.

### Step 3: Auto-Groom (if needed)

Check if a grooming manifest exists for this item:

1. Search `.claude/grooming-reports/` for files containing the item title
2. Search conversation context for a recent `/groom-backlog-item` output matching this item

If no grooming context exists:

```text
Skill(skill="groom-backlog-item", args="{item title}")
```

Capture the grooming output (context manifest with Related Research, Supporting Skills, Related Agents, Prior Work, Dependencies, Blockers, Suggested First Steps).

### Step 4: Compose Feature Request

Build the feature request string that will become `$ARGUMENTS` for `add-new-feature`:

```text
## Backlog Item: {title}

**Source**: {source}
**Priority**: {priority section - P0/P1/P2/Ideas}
**Added**: {added date}

### Description

{description text}

### Research Questions

{research_first text, or "None" if absent}

### Suggested Location

{suggested_location text, or "To be determined during architecture phase" if absent}

### Grooming Context

{full context manifest from Step 3, if available}
```

### Step 5: Invoke SAM Planning

```text
Skill(skill="python3-development:add-new-feature", args="{composed feature request}")
```

This runs the full SAM workflow: discovery, codebase analysis, architecture spec, task decomposition, validation, context manifest.

### Step 6: Update Backlog with Plan Reference

After `add-new-feature` completes, identify the task file it created by searching:

```text
Glob(pattern="plan/tasks-*-{slug}*")
```

Where `{slug}` is derived from the item title (lowercased, spaces to hyphens).

Add a `**Plan**:` field to the backlog item in `.claude/BACKLOG.md`:

```text
### SAM: {item title}

**Source**: {source}
**Added**: {added date}
**Plan**: plan/tasks-{N}-{slug}.md
**Description**: {description}
```

Update `last-updated:` in the BACKLOG.md YAML frontmatter to today's date.

### Step 7: Report Next Steps

```text
Backlog item "{title}" is now planned.

- Plan file: plan/tasks-{N}-{slug}.md (or plan/tasks-{N}-{slug}/ directory)
- To execute: /python3-development:implement-feature {slug}
- To check status: /python3-development:implementation-manager status . {slug}
```

## Error Handling

- Item not found: list available items from BACKLOG.md with their priority sections
- Multiple matches: present numbered list, ask user to choose
- Grooming fails: proceed without grooming context, note the gap in the feature request
- add-new-feature fails: report the failure, do not update BACKLOG.md
- Plan file not found after add-new-feature: search `plan/` directory broadly, ask user to confirm the path

## Example Session

```text
> /work-backlog-item Error Recovery

Found: "SAM: Error Recovery / Rollback Procedures" (P1)
No grooming manifest found. Running /groom-backlog-item first...

[grooming output]

Composing feature request from backlog item + grooming context...
Invoking /python3-development:add-new-feature...

[SAM phases run]

Updated BACKLOG.md with Plan: plan/tasks-2-error-recovery.md

Next steps:
- To execute: /python3-development:implement-feature error-recovery
- To check status: /python3-development:implementation-manager status . error-recovery
```
