---
description: Bridges BACKLOG.md to the SAM planning pipeline — use when you want to pick a backlog item and move it into a SAM plan. No args shows interactive backlog browser with grooming status; with args finds item by title substring, auto-grooms if needed, runs RT-ICA to BLOCK on missing inputs before SAM planning, invokes add-new-feature, then updates backlog with plan reference. STOPS if item already has a Plan field or RT-ICA returns BLOCKED.
argument-hint: '[item-title-substring]'
user-invocable: true
---

# Work Backlog Item

Bridge a `.claude/BACKLOG.md` item into the SAM planning pipeline via `/python3-development:add-new-feature`, then record the resulting plan file back in BACKLOG.md.

When invoked with no arguments, shows an interactive backlog browser. When invoked with a title substring, proceeds directly to the planning workflow.

## Arguments

`$ARGUMENTS` is an optional case-insensitive title substring matching an H3 heading in BACKLOG.md.

```text
/work-backlog-item                    # interactive browser
/work-backlog-item Error Recovery     # direct match
/work-backlog-item regex false positive
/work-backlog-item validate-task-file
```

## Workflow

### Step 0: Interactive Browser (no arguments only)

**Trigger:** `$ARGUMENTS` is empty.

<step0_procedure>

1. Read `.claude/BACKLOG.md`. Parse all H3 headings (`### ...`) from P0, P1, P2, and Ideas sections. Record each item's priority section and title.

2. For each item, determine grooming status:
   - **Has plan** — item has a `**Plan**:` field in BACKLOG.md
   - **Groomed** — `.claude/grooming-reports/` contains a file whose content references the item title (search with Grep)
   - **Ungroomed** — neither condition above is true

3. Present a numbered list. Use these status indicators in user-visible output only:

   ```text
   Backlog Items:

   P0
     1. ✅ SAM: Error Recovery / Rollback Procedures
     2. 🔍 SAM: Regex False Positive Suppression

   P1
     3. 📋 SAM: Validate Task File Schema
     4. 📋 SAM: Implement Feature Dry-Run Mode

   P2
     5. 🔍 SAM: Context Window Budget Tracking

   Ideas
     6. 📋 SAM: Multi-Repo Support

   Status: ✅ = planned  🔍 = groomed  📋 = not yet groomed

   Options:
     [number]   — Select item to work on
     G [number] — Groom a specific item
     G all      — Groom all ungroomed items
     D [number] — Show full details for an item
   ```

4. Use `AskUserQuestion` to ask: "Which item would you like to work on next?"

5. Handle the response:
   - `[number]` — set `$ARGUMENTS` to that item's title and proceed to Step 1
   - `G [number]` — invoke `Skill(skill="groom-backlog-item", args="{item title}")` then re-display the list
   - `G all` — invoke `Skill(skill="groom-backlog-item", args="all")` then re-display the list
   - `D [number]` — display the full item description, research_first field, and grooming manifest if it exists in `.claude/grooming-reports/`, then re-display the list

</step0_procedure>

### Step 1: Find the Backlog Item

Read `.claude/BACKLOG.md`. Search H3 headings (`### ...`) for case-insensitive match against `$ARGUMENTS`.

- Zero matches: report "No backlog item found matching: {$ARGUMENTS}" and stop.
- Multiple matches: list all matches and ask the user to pick one.

Record the priority section (P0, P1, P2, Ideas) the item belongs to.

### Step 2: Extract Item Fields

From the matched item, extract these fields (all are `**bold-key**: value` format on separate lines):

- `title` — the H3 heading text (required)
- `source` — `**Source**:` value (optional)
- `added` — `**Added**:` value (optional)
- `description` — `**Description**:` value (required)
- `research_first` — `**Research first**:` value (optional)
- `suggested_location` — `**Suggested location**:` value (optional)
- `plan` — `**Plan**:` value (optional)

If the item already has a `**Plan**:` field, report:

```text
This item already has a plan at {path}. Use /python3-development:implement-feature {path} to execute it.
```

Then stop.

### Step 3: Auto-Groom (if needed)

<groom_check>

1. Search `.claude/grooming-reports/` for files whose content references the item title (use Grep)
2. Search conversation context for a recent `groom-backlog-item` output matching this item.

If no grooming report exists:

```text
Skill(command: "groom-backlog-item", args: "{item title}")
```

Capture the grooming output — context manifest with Related Research, Supporting Skills, Related Agents, Prior Work, Dependencies, Blockers, Suggested First Steps.
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

After `add-new-feature` completes, identify the task file it created:

```text
Glob(pattern="plan/tasks-*-{slug}*")
```

Where `{slug}` is the item title lowercased with spaces replaced by hyphens.

Add a `**Plan**:` field to the backlog item in `.claude/BACKLOG.md`:

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

- Item not found: list available items from BACKLOG.md with their priority sections
- Multiple matches: present numbered list, ask user to choose
- Grooming fails: proceed without grooming context, note the gap in the feature request
- RT-ICA returns BLOCKED: present missing inputs, wait for user, do not invoke `add-new-feature`
- `add-new-feature` fails: report the failure, do not update BACKLOG.md
- Plan file not found after planning: search `plan/` directory broadly, ask user to confirm the path
- Grooming reports directory does not exist: treat all items as ungroomed

## Example Session (with argument)

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
