---
name: work-backlog-item
description: "Use when working, planning, or closing a backlog item. Bridges BACKLOG.md to the SAM planning pipeline with optional GitHub Issue/Project/Milestone tracking. No args: interactive browser. With '#N': load item directly from GitHub Issue #N (labels and milestone are canonical status). With title substring: auto-grooming, RT-ICA gate, GitHub issue sync, SAM planning, plan reference recorded. '--auto {title}': fully autonomous mode — no AskUserQuestion calls, derives missing data from research files, logs all decisions, skips interactive GitHub prompts; suitable for agent use without human in the loop. 'close {title}': verifies plan checklist 100% complete, closes GitHub issue, marks DONE. 'resolve {title}': marks item no longer applicable with reason. 'setup-github': initializes labels, creates project and first milestone. STOPS if item has existing Plan field or RT-ICA returns BLOCKED."
argument-hint: '[#N | --auto {title} | item-title-substring | close {title} | resolve {title} | setup-github]'
user-invocable: true
---
# Work Backlog Item

Bridge a backlog item into the SAM planning pipeline via `/python3-development:add-new-feature`.
Primary source of truth is **GitHub Issues** (labels + milestone = canonical status); `.claude/BACKLOG.md` is the local scratchpad and is kept in sync.

When invoked with no arguments, shows an interactive browser. When invoked with `#N` or a title substring, proceeds directly to the planning workflow.

## Arguments

`$0` selects the operating mode; remaining positional args (`$1`, `$2`, ...) form the title or parameter:

| `$0` value | Remaining args | Mode |
|---|---|---|
| (empty) | — | Interactive browser |
| `#N` | — | Issue-first: load item from GitHub Issue #N |
| `--auto` | `$1`+ = title (or empty → auto-select first open P0/P1 item) | Autonomous — no `AskUserQuestion` calls |
| `close` | `$1`+ = title or `#N` | Verify and close a completed item |
| `resolve` | `$1`+ = title or `#N` | Mark no longer applicable (reason required) |
| `setup-github` | — | Initialize labels, project, first milestone |
| (any other) | — | `$ARGUMENTS` treated as title substring → planning |

```text
/work-backlog-item                                    # interactive browser
/work-backlog-item #42                               # issue-first → planning
/work-backlog-item Error Recovery                    # direct match → planning
/work-backlog-item --auto                            # autonomous → auto-select first open P0/P1
/work-backlog-item --auto vercel skills npm package  # autonomous → planning
/work-backlog-item close Error Recovery              # verify and close by title
/work-backlog-item close #42                         # verify and close by issue number
/work-backlog-item resolve commitlint                # mark no longer applicable
/work-backlog-item resolve #17                       # mark no longer applicable by issue
```

### --auto mode rules

When `$0` is `--auto`, the following substitutions apply at every interactive decision point:

| Normal behaviour | `--auto` substitution |
|---|---|
| No title given (`$1` is empty) | Scan BACKLOG.md P0 then P1 sections for the first item that has no `**Completed**:` line, no `**Status**: DONE`, and no ~~strikethrough~~ heading. Log `[AUTO] No title — auto-selected: {title}`, use that item's title as the working title and proceed. If no open P0/P1 item exists, log `[AUTO] STOP — no open P0/P1 items found` and stop. |
| Step 1b: issue not found | Log `[AUTO] STOP — Issue #N not found`, stop |
| Step 1: zero matches → ask user to create | Auto-invoke `create-backlog-item --auto {title}`, log `[AUTO] No item found — invoking create-backlog-item --auto` |
| Step 1: multiple matches → ask user to pick | Log `[AUTO] Multiple matches — picking first: {title}`, proceed with first match |
| Step 2.5: offer GitHub issue for P0/P1 | Log `[AUTO] Skipping GitHub issue offer`, continue without issue |
| Step 2.5: ask milestone assignment | Log `[AUTO] Skipping milestone assignment`, skip |
| RT-ICA BLOCKED | Log `[AUTO] STOP — RT-ICA BLOCKED: {missing inputs}`, stop (cannot resolve without human) |
| Any other `AskUserQuestion` | Log `[AUTO] Decision: {chosen option} — reason: {evidence}`, proceed with logged choice |

`--auto` does NOT change the behaviour of Steps 3–8 (grooming, RT-ICA evaluation, SAM planning, BACKLOG.md write) — those are already agent-executable without human input.

## Workflow

### Routing (evaluated first, before any step)

Dispatch based on `$0` (the first argument word) before executing any step:

| `$0` value | Title source | Route |
|---|---|---|
| (empty) | — | Step 0 — interactive browser |
| `#N` (starts with `#`) | issue number | Step 1b — Issue-first path |
| `--auto` | `$1`+ joined (empty → auto-select first open P0/P1) | AUTO_MODE=true → Step 1 |
| `close` | `$1`+ joined (title or `#N`) | Step 9 (close path) |
| `resolve` | `$1`+ joined (title or `#N`) | Step 9 (resolve path) |
| `setup-github` | — | setup-github command |
| (any other) | `$ARGUMENTS` | Title substring → Step 1 (interactive mode) |

**AUTO_MODE** — when set, all `AskUserQuestion` calls are replaced with evidence-derived decisions. See the `--auto mode rules` table in the Arguments section for each substitution.

### Step 0: Interactive Browser (no arguments only)

**Trigger:** `$0` is empty (no arguments passed).

<step0_procedure>

1. Read `.claude/BACKLOG.md`. Parse all H3 headings (`### ...`) from P0, P1, P2, and Ideas sections. Record each item's priority section and title.

2. For each item, determine status using two sources (prefer GitHub when an `**Issue**: #N` field is present):
   - **GitHub status** (preferred when `**Issue**: #N` exists) — fetch `status:*` label and milestone title from:
     ```bash
     gh issue view {issue_number} -R Jamie-BitFlight/claude_skills --json number,state,labels,milestone \
       --jq '[.state, (.labels|map(.name)|join(",")), (.milestone.title // "—")] | @tsv'
     ```
   - **Local status** (fallback when no `**Issue**: #N`) — check for `**Plan**:` field (planned) or grooming report in `.claude/grooming-reports/` (groomed)

3. Present a numbered list. Use these status indicators in user-visible output only:

   ```text
   Backlog Items:

   P0
     1. ✅ SAM: Error Recovery / Rollback Procedures       [#12  status:in-progress  v1.0]
     2. 🔍 SAM: Regex False Positive Suppression           [#14  status:needs-grooming  v1.0]

   P1
     3. 📋 SAM: Validate Task File Schema                  [no issue]
     4. 📋 SAM: Implement Feature Dry-Run Mode             [no issue]

   P2
     5. 🔍 SAM: Context Window Budget Tracking             [#18  status:needs-grooming  —]

   Ideas
     6. 📋 SAM: Multi-Repo Support                         [no issue]

   Status: ✅ = planned/in-progress  🔍 = groomed/needs-grooming  📋 = not yet groomed

   Options:
     [number]   — Select item to work on
     G [number] — Groom a specific item
     G all      — Groom all ungroomed items
     D [number] — Show full details for an item
   ```

4. Use `AskUserQuestion` to ask: "Which item would you like to work on next?"

5. Handle the response:
   - `[number]` — use that item's title as the working title and proceed to Step 1
   - `G [number]` — invoke `Skill(skill="groom-backlog-item", args="{item title}")` then re-display the list
   - `G all` — invoke `Skill(skill="groom-backlog-item", args="all")` then re-display the list
   - `D [number]` — display the full item description, research_first field, and grooming manifest if it exists in `.claude/grooming-reports/`, then re-display the list
   - `C [number]` — proceed to Step 9 (close path) with that item's title
   - `R [number]` — proceed to Step 9 (resolve path) with that item's title

</step0_procedure>

**Routing:** If `$0` is `close` or `resolve`, extract `$1`+ as the title and jump directly to Step 9.

### Step 1b: Issue-First Path (`#N`)

**Trigger:** `$0` matches `#[0-9]+`.

<issue_first_procedure>

Extract the issue number `{issue_number}` from `$0`. Fetch the issue:

```bash
gh issue view {issue_number} -R Jamie-BitFlight/claude_skills \
  --json number,title,state,body,labels,milestone
```

If the issue does not exist, report and stop.
If the issue is already closed, warn: "Issue #{issue_number} is closed. Use `close` or `resolve` if needed." and stop.

From the issue response build the working item:

| Field | Source |
|---|---|
| `title` | issue `title` |
| `description` | issue `body` (full text) |
| `source` | `"GitHub Issue #N"` |
| `priority` | `priority:*` label → P0 / P1 / P2 / Ideas |
| `status` | `status:*` label (canonical — do not read BACKLOG.md for status) |
| `milestone` | issue `milestone.title` |
| `plan` | search `body` for `**Plan**:` line |

Then try to find a matching item in `.claude/BACKLOG.md` by issue number (`**Issue**: #N`) or title. If found, use it to supplement any missing fields (e.g. `research_first`, `suggested_location`). If not found, continue without a BACKLOG.md record — the GitHub Issue is sufficient.

Skip to Step 3 with the assembled item.

</issue_first_procedure>

### Step 1: Find the Backlog Item

Read `.claude/BACKLOG.md`. Search H3 headings (`### ...`) for case-insensitive match against the title. Title = `$1`+ joined (args after the mode flag `$0`). In interactive mode, title = full `$ARGUMENTS`.

**AUTO_MODE with no title (`$1` is empty):** apply the "No title given" substitution from the `--auto mode rules` table — scan P0 then P1 sections for the first open item, log and use its title. Skip items whose H3 heading contains ~~strikethrough~~, or whose body contains a `**Completed**:` line or `**Status**: DONE`.

- **Zero matches (interactive mode):** report "No backlog item found matching: {title}" and offer to create one via `/create-backlog-item`.
- **Zero matches (AUTO_MODE):** log `[AUTO] No item found — invoking create-backlog-item --auto {title}`, invoke `Skill(command: "create-backlog-item", args: "--auto {title}")`, then re-run Step 1.
- **Multiple matches (interactive mode):** list all matches and ask the user to pick one.
- **Multiple matches (AUTO_MODE):** log `[AUTO] Multiple matches — picking first: {title}`, proceed with first match.

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

Build the feature request string for `add-new-feature`:

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

If the item has `**Issue**: #N`, record it in the plan file header comment and include `Fixes #N` in any commit message produced during implementation.

Update `last-updated:` in the BACKLOG.md YAML frontmatter to today's date.

### Step 8: Report Next Steps

```text
Backlog item "{title}" is now planned.

- Plan file: plan/tasks-{N}-{slug}.md (or plan/tasks-{N}-{slug}/ directory)
- To execute:      /python3-development:implement-feature {slug}
- To check status: /python3-development:implementation-manager status . {slug}
- To close when done: /work-backlog-item close {slug}
```

### Step 9: Verify and Close

**Trigger:** `$0` is `close` or `resolve`.

<step9_procedure>

Extract the operation from `$0` and the argument from `$1`+:

- `$0` = `close`: `$1`+ = title or `#N` → verify implementation and mark COMPLETED
- `$0` = `resolve`: `$1`+ = title or `#N` → mark no longer applicable (no verification required)

#### 9a: Find Item

If `$1` starts with `#` (e.g., `close #42`), treat it as an issue number:

```bash
gh issue view {issue_number} -R Jamie-BitFlight/claude_skills \
  --json number,title,state,body,labels
```

- If the issue is not found, report and stop.
- Extract `title` from the issue response and use it as the working title.

Otherwise, read `.claude/BACKLOG.md` and search H3 headings for case-insensitive match against `{title}`.

- Zero matches: report "No backlog item found matching: {title}" and stop.
- Multiple matches: list all matches and ask user to pick one.
- Item already in `## Completed` section: report "Item already closed on {Completed date}" and stop.

#### 9b: Resolve path (skip verification)

If operation is `resolve`:

1. Use `AskUserQuestion` to ask: "Why is this item no longer applicable?" (free text)
2. Update the matched item in `.claude/BACKLOG.md`:

```text
### {original title}

**Source**: {source}
**Added**: {added}
**Resolved**: {YYYY-MM-DD}
**Status**: RESOLVED — {user-provided reason}
**Description**: {description}
```

3. Update `last-updated:` in BACKLOG.md YAML frontmatter to today's date.
4. Report:

```text
Backlog item "{title}" resolved.
Reason: {reason}
```

Then stop.

#### 9c: Close path — checklist verification

If operation is `close`:

1. Extract `**Plan**:` field from the matched item. If absent:

```text
No plan file recorded for "{title}". Cannot verify checklist.
Either run /work-backlog-item {title} first to create a plan,
or use /work-backlog-item resolve {title} if no plan was needed.
```

Then stop.

2. Read the plan file. Count:
   - `total_tasks` — lines matching `- \[ \]` or `- \[x\]`
   - `checked_tasks` — lines matching `- \[x\]`

3. If `checked_tasks < total_tasks`:

```text
Checklist incomplete: {checked_tasks}/{total_tasks} tasks done.

Remaining:
{list of unchecked task lines}

Complete all tasks before closing this item.
```

Then stop.

#### 9d: Close path — acceptance criteria verification

4. Spawn a verification agent:

```text
Task(
  subagent_type: "general-purpose",
  prompt: "You are verifying whether a completed backlog item genuinely satisfies its stated goal.

Backlog item title: {title}
Description and acceptance criteria:
{description text from BACKLOG.md}

Plan file: {plan file path}
Plan checklist: {checked_tasks}/{total_tasks} — 100% complete.

Your task:
1. Read the plan file to understand what was implemented.
2. Search git log for commits referencing this item (use: git log --oneline -20).
3. Read 2-3 key changed files to verify the implementation exists.
4. Assess: Does the implementation satisfy the stated goal? Is the product better for it?

Return:
- PASS or FAIL
- One sentence of evidence (file:line or commit SHA)
- Any gaps you found (if FAIL)"
)
```

5. Collect agent verdict:
   - **PASS**: proceed to 9e
   - **FAIL**: report gaps, do not close:

```text
Verification FAILED for "{title}".

Gaps found:
{agent findings}

Address these gaps before closing.
```

Then stop.

#### 9e: Write closing record

6. Update the matched item in `.claude/BACKLOG.md`:

```text
### {original title}

**Source**: {source}
**Added**: {added}
**Completed**: {YYYY-MM-DD}
**Status**: DONE — verified by checklist ({checked}/{total}) + acceptance criteria check
**Plan**: {plan file path}
**Description**: {description}
```

7. Update `last-updated:` and `last-completed:` in BACKLOG.md YAML frontmatter to today's date.

8. Report:

```text
Backlog item "{title}" closed.

- Checklist: {checked}/{total} tasks complete
- Acceptance criteria: PASS
- Status written to BACKLOG.md
```

</step9_procedure>

## GitHub Integration

BACKLOG.md is the local development scratchpad. GitHub Issues are the published, tracked view. They are linked — when you work a P0/P1 backlog item, a GitHub Issue can be created and kept in sync.

```text
BACKLOG.md          →  GitHub Issue
  Priority section  →  priority:* label
  Description       →  Issue body (story format)
  Status            →  status:* label
  Plan file         →  Issue body Notes
  **Issue**: #N     ←  written back after creation
  Completed         →  Issue closed
```

Full step-by-step commands and example sessions: [github-integration.md](./references/github-integration.md)

### Commit Messages and PR Body — `Fixes #N`

When implementing a backlog item that has a linked GitHub Issue (`**Issue**: #N`), **every commit and the PR description MUST reference the issue** so GitHub auto-closes it on merge:

- **Commit messages**: end the subject line with `(fixes #N)` or add a `Fixes #N` footer in the body:
  ```text
  fix(perl-development): replace shell injection in run_command (fixes #80)
  ```
  or
  ```text
  fix(perl-development): replace shell injection in run_command

  Fixes #80
  ```
- **PR description**: include at least one `Fixes #N` line in the PR body so GitHub links and auto-closes the issue when the PR is merged.
- When a PR fixes **multiple issues**, list each on a separate line:
  ```text
  Fixes #80
  Fixes #83
  ```
- Use `Fixes` (not `Closes`) for bugs and security issues; either works for feature work.

This convention ensures issues are automatically closed on merge without manual intervention.

### Step 2.5: GitHub Issue Sync

After Step 2, check for `**Issue**: #N` field in the matched item.

- Found: verify issue state with `gh issue view {issue_number} -R Jamie-BitFlight/claude_skills --json number,title,state,labels`
- Not found + P0/P1: offer to create a GitHub Issue (proceed to Step 2.5a)
- Not found + P2/Ideas: skip silently

### Step 2.5a: Create GitHub Issue

Build story-format body (Story / Description / Acceptance Criteria / Context sections). Use the Python script to create the issue with the correct labels:

```bash
uv run .claude/skills/gh/scripts/github_project_setup.py issue create \
  --repo Jamie-BitFlight/claude_skills \
  --title "{conventional-commits-type}: {item title}" \
  --body "{story body}" \
  --priority-label "priority:{p0|p1|p2|idea}" \
  --type-label "type:{feature|bug|refactor|docs|chore}" \
  --milestone {milestone_number}
```

Capture the issue number from output and write `**Issue**: #N` back to BACKLOG.md.

### Step 2.7: Set In-Progress Label

If the item has `**Issue**: #N`, transition the issue to in-progress:

```bash
uv run .claude/skills/gh/scripts/github_project_setup.py milestone start \
  --number {milestone_number} --repo Jamie-BitFlight/claude_skills
```

If no milestone is assigned, transition only the single issue's label directly:

```bash
gh issue edit {issue_number} -R Jamie-BitFlight/claude_skills \
  --add-label "status:in-progress" \
  --remove-label "status:needs-grooming"
```

### Step 9 Extension: Close GitHub Issue

After writing the closing record (Step 9e), determine the issue number:

- If invoked as `close #N` or `resolve #N`: use the issue number from `$1`.
- Otherwise, check the matched BACKLOG.md item for an `**Issue**: #N` field.

If an issue number is found:

```bash
gh issue close {issue_number} -R Jamie-BitFlight/claude_skills \
  --comment "Completed. Checklist {checked}/{total} — PASS. Plan: {plan file path}"
```

If no issue number is available, skip silently.

### setup-github Command

**Trigger:** `$0` is `setup-github`. Initializes label taxonomy, first milestone, and GitHub Project.

```bash
uv run .claude/skills/gh/scripts/github_project_setup.py labels --repo Jamie-BitFlight/claude_skills
gh api repos/Jamie-BitFlight/claude_skills/milestones -X POST \
  -f title="v1.0 — Skills Foundation" -f due_on="2026-03-31T00:00:00Z"
gh project create --owner Jamie-BitFlight --title "claude_skills Backlog"
```

Full setup steps and expected output: [github-integration.md](./references/github-integration.md)

## Error Handling

- `#N` not found: report and list open issues with `gh issue list -R Jamie-BitFlight/claude_skills --state open`
- `#N` already closed: warn and stop; offer `close` or `resolve` if needed
- `close #N` / `resolve #N` — issue not found: report and stop
- Item not found: list available items from BACKLOG.md with their priority sections
- Multiple matches: present numbered list, ask user to choose
- Grooming fails: proceed without grooming context, note the gap in the feature request
- RT-ICA returns BLOCKED: present missing inputs, wait for user, do not invoke `add-new-feature`
- `add-new-feature` fails: report the failure, do not update BACKLOG.md
- Plan file not found after planning: search `plan/` directory broadly, ask user to confirm the path
- Grooming reports directory does not exist: treat all items as ungroomed
- `close` with no `**Plan**:` field: report and offer `resolve` as alternative
- `close` with incomplete checklist: list remaining tasks, do not close
- `close` with verification FAIL: report gaps, do not close
- `close` on already-completed item: report closed date, do not re-close
- `resolve` with no reason provided: block until user provides reason (reason is required evidence)
- GitHub issue creation fails: report error, continue with BACKLOG.md-only workflow; do not block SAM planning
- `gh` not installed: run `uv run .claude/skills/gh/scripts/setup_gh.py` first
- Label not found during issue create: `github_project_setup.py` creates it automatically
- Milestone not found: skip milestone assignment; do not fail

## Example Sessions

### Issue-first planning (`#N`)

```text
> /work-backlog-item #131

Loading GitHub Issue #131...
  Title:     plugin-validator: UX and coverage gaps
  Labels:    priority:p1, status:needs-grooming
  Milestone: v1.1 — Quality Gates
  State:     open

Matched BACKLOG.md item for additional context: ✓
No grooming manifest found. Running groom-backlog-item first...

[grooming output]

RT-ICA: APPROVED — all conditions available.
Setting status:in-progress on issue #131...
  ✓ status:needs-grooming → status:in-progress

Composing feature request...
Invoking /python3-development:add-new-feature...

[SAM phases run]

Updated BACKLOG.md with Plan: plan/tasks-2-validator-ux-coverage.md

Next steps:
- To execute:      /python3-development:implement-feature validator-ux-coverage
- To close when done: /work-backlog-item close plugin-validator UX and coverage gaps
```

### Planning (with title substring)

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
- To close when done: /work-backlog-item close error-recovery
```

### Closing a completed item (by title)

```text
> /work-backlog-item close validator UX

Found: "plugin-validator UX and coverage gaps" (P1)
Plan: plan/tasks-2-validator-ux-coverage.md
Checklist: 12/12 tasks complete

Spawning acceptance criteria verification agent...

Verdict: PASS
Evidence: Sub-issues 1-4 implemented in plugins/plugin-creator/scripts/plugin_validator.py
          commit 4a2f1b3 — "fix(validator): report unique files, add hook validation"

Backlog item "plugin-validator UX and coverage gaps" closed.
- Checklist: 12/12 tasks complete
- Acceptance criteria: PASS
- Status written to BACKLOG.md
```

### Closing a completed item (by issue number)

```text
> /work-backlog-item close #131

Fetching GitHub Issue #131...
  Title: plugin-validator UX and coverage gaps
  State: open

Found BACKLOG.md match: "plugin-validator UX and coverage gaps" (P1)
Plan: plan/tasks-2-validator-ux-coverage.md
Checklist: 12/12 tasks complete

Spawning acceptance criteria verification agent...

Verdict: PASS
Evidence: commit 4a2f1b3 — "fix(validator): report unique files, add hook validation"

Backlog item "plugin-validator UX and coverage gaps" closed.
- Checklist: 12/12 tasks complete
- Acceptance criteria: PASS
- Status written to BACKLOG.md
- GitHub Issue #131 closed
```

### Resolving a no-longer-applicable item

```text
> /work-backlog-item resolve commitlint verify last flag

Found: "commitlint: Verify --last flag and exit codes against primary sources" (P1)
Why is this item no longer applicable?
> REFUTED by fact-check: --last flag is verified in commitlint source cli.ts. No fix needed.

Backlog item resolved.
  Resolved: 2026-02-21
  Status: RESOLVED — REFUTED by fact-check: --last flag verified against commitlint
          source cli.ts and official docs. No fix needed.
```

GitHub-specific example sessions (issue creation flow and setup-github): [github-integration.md](./references/github-integration.md)

---

## Validation Plan

See [validation-plan.md](./references/validation-plan.md) for V1–V6 verification commands and the full integration test sequence.
