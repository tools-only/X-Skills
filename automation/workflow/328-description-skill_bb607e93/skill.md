---
description: "Invoke when the user starts work on a backlog item — TRIGGER item title substring provided; reads .claude/BACKLOG.md, auto-grooms if no manifest exists, runs RT-ICA to BLOCK on missing inputs before SAM planning, invokes add-new-feature, then writes Plan reference back to BACKLOG.md. STOPS if item already has a Plan field or RT-ICA returns BLOCKED."
argument-hint: <item-title-substring>
user-invocable: true
---

# Work Backlog Item

Bridge a `.claude/BACKLOG.md` item into the SAM planning pipeline via `/python3-development:add-new-feature`, then record the resulting plan file back in BACKLOG.md.

## Arguments

`$ARGUMENTS` — case-insensitive substring matching an H3 heading in BACKLOG.md.

```text
/work-backlog-item Error Recovery
/work-backlog-item regex false positive
/work-backlog-item validate-task-file
```

## Workflow

### Step 1: Find the Backlog Item

Read `.claude/BACKLOG.md`. Search H3 headings (`### ...`) for case-insensitive match against `$ARGUMENTS`.

- Zero matches: report `No backlog item found matching: {$ARGUMENTS}` and stop.
- Multiple matches: list all matches and ask the user to pick one.

Record the priority section (P0, P1, P2, Ideas) the item belongs to.

### Step 2: Extract Item Fields

Extract these fields (all use **bold-key**: value format on separate lines):

| Field              | Key in BACKLOG.md         | Required |
|--------------------|---------------------------|----------|
| title              | H3 heading text           | Yes      |
| source             | `**Source**:`             | No       |
| added              | `**Added**:`              | No       |
| description        | `**Description**:`        | Yes      |
| research_first     | `**Research first**:`     | No       |
| suggested_location | `**Suggested location**:` | No       |
| plan               | `**Plan**:`               | No       |

If the item already has a `**Plan**:` field, report:

```text
This item already has a plan at {path}. Use /python3-development:implement-feature {path} to execute it.
```

Then stop.

### Step 3: Auto-Groom (if needed)

<groom_check>
1. Search `.claude/grooming-reports/` for files containing the item title.
2. Search conversation context for a recent `groom-backlog-item` output matching this item.

If no grooming context exists, invoke:

```text
Skill(command: "groom-backlog-item", args: "{item title}")
```

Capture the grooming output: context manifest containing Related Research, Supporting Skills, Related Agents, Prior Work, Dependencies, Blockers, Suggested First Steps.
</groom_check>

### Step 4: RT-ICA Checkpoint

<rtica_gate>
Before composing the feature request, verify the grooming context manifest contains an RT-ICA summary. If absent, perform RT-ICA now:

1. **Goal statement** — What completing this item achieves.
2. **Reverse prerequisites** — Conditions required for success (enumerate each).
3. **Availability check** — For each condition: AVAILABLE / DERIVABLE / MISSING.
4. **Decision** — APPROVED or BLOCKED.

**If BLOCKED:**

Present a structured summary to the user:

```text
RT-ICA: BLOCKED

Missing inputs that prevent SAM planning:
- {missing condition 1}
- {missing condition 2}

Provide these inputs before proceeding. SAM planning will not be invoked with known gaps.
```

Wait for user response. Do not invoke Step 5 until BLOCKED is resolved.

**If APPROVED:**

Proceed to Step 5. Carry DERIVABLE items forward as "Assumptions to confirm" in the RT-ICA section of the feature request.
</rtica_gate>

### Step 5: Compose Feature Request

Build the `$ARGUMENTS` string for `add-new-feature`:

```text
## Backlog Item: {title}

**Source**: {source}
**Priority**: {priority section — P0/P1/P2/Ideas}
**Added**: {added date}

### Description

{description text}

### Research Questions

{research_first text, or "None" if absent}

### Suggested Location

{suggested_location text, or "To be determined during architecture phase" if absent}

### RT-ICA Assessment

**Decision**: APPROVED
**Goal**: {goal statement}
**Verified conditions**: {list of AVAILABLE items}
**Assumptions to confirm**: {list of DERIVABLE items, or "None"}

### Grooming Context

{full context manifest from Step 3, if available}
```

### Step 6: Invoke SAM Planning

```text
Skill(command: "python3-development:add-new-feature", args: "{composed feature request}")
```

This runs the full SAM workflow: discovery, codebase analysis, architecture spec, task decomposition, validation, context manifest.

### Step 7: Update Backlog with Plan Reference

After `add-new-feature` completes, derive the task file slug:

- slug = item title lowercased, spaces replaced with hyphens

Search for the generated plan file:

```text
Glob(pattern="plan/tasks-*-{slug}*")
```

Add a `**Plan**:` field to the matched item in `.claude/BACKLOG.md`:

```text
### SAM: {item title}

**Source**: {source}
**Added**: {added date}
**Plan**: plan/tasks-{N}-{slug}.md
**Description**: {description}
```

Update `last-updated:` in the BACKLOG.md YAML frontmatter to today's date.

### Step 8: Report Next Steps

```text
Backlog item "{title}" is now planned.

- Plan file: plan/tasks-{N}-{slug}.md (or plan/tasks-{N}-{slug}/ directory)
- To execute:      /python3-development:implement-feature {slug}
- To check status: /python3-development:implementation-manager status . {slug}
```

## Error Handling

| Condition                          | Action                                                                  |
|------------------------------------|-------------------------------------------------------------------------|
| Item not found                     | List available items from BACKLOG.md with priority sections             |
| Multiple matches                   | Present numbered list, ask user to choose                               |
| Grooming fails                     | Proceed without grooming context; note the gap in the feature request   |
| RT-ICA returns BLOCKED             | Present missing inputs; wait for user; do not invoke `add-new-feature` |
| `add-new-feature` fails            | Report the failure; do not update BACKLOG.md                            |
| Plan file not found after planning | Search `plan/` directory broadly; ask user to confirm the path          |

## Example Session

```text
> /work-backlog-item Error Recovery

Found: "SAM: Error Recovery / Rollback Procedures" (P1)
No grooming manifest found. Running groom-backlog-item first...

[grooming output]

RT-ICA: APPROVED — all conditions available.
Composing feature request...
Invoking /python3-development:add-new-feature...

[SAM phases run]

Updated BACKLOG.md with Plan: plan/tasks-2-error-recovery.md

Next steps:
- To execute:      /python3-development:implement-feature error-recovery
- To check status: /python3-development:implementation-manager status . error-recovery
```
