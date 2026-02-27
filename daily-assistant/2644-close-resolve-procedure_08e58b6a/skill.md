# Close / Resolve Procedure

Full Step 9 workflow. Trigger: `$0` is `close` or `resolve`.

Extract the operation from `$0` and the argument from `$1`+:

- `$0` = `close`: `$1`+ = title or `#N` → verify implementation and mark COMPLETED
- `$0` = `resolve`: `$1`+ = title or `#N` → mark no longer applicable (no verification required)

## 9a: Find Item

If `$1` starts with `#` (e.g., `close #42`), treat it as an issue number:

```bash
gh issue view {issue_number} -R Jamie-BitFlight/claude_skills \
  --json number,title,state,body,labels
```

- If the issue is not found, report and stop.
- Extract `title` from the issue response and use it as the working title.

Otherwise, scan `.claude/backlog/` per-item files and search item titles for case-insensitive match against `{title}`.

- Zero matches: report "No backlog item found matching: {title}" and stop.
- Multiple matches: list all matches and ask user to pick one.
- Item already in `## Completed` section: report "Item already closed on {Completed date}" and stop.

## 9b: Resolve path (skip verification)

If operation is `resolve`:

1. Use `AskUserQuestion` to ask: "Why is this item no longer applicable?" (free text)
2. Invoke the backlog script:

```bash
uv run .claude/skills/backlog/scripts/backlog.py resolve "{title or #N}" --reason "{reason}" -R Jamie-BitFlight/claude_skills
```

3. Report the script output to the user.

Then stop.

## 9c: Close path — checklist verification

If operation is `close`:

1. Extract `**Plan**:` field from the matched item. If absent:

<eg>
No plan file recorded for "{title}". Cannot verify checklist.
Either run /work-backlog-item {title} first to create a plan,
or use /work-backlog-item resolve {title} if no plan was needed.
</eg>

Then stop.

2. Read the plan file. Count:
   - `total_tasks` — lines matching `- \[ \]` or `- \[x\]`
   - `checked_tasks` — lines matching `- \[x\]`

3. If `checked_tasks < total_tasks`:

<eg>
Checklist incomplete: {checked_tasks}/{total_tasks} tasks done.

Remaining:
{list of unchecked task lines}

Complete all tasks before closing this item.
</eg>

Then stop.

## 9d: Close path — typed acceptance-criteria verification

4. Extract `**Acceptance Criteria**:` from the per-item file (`.claude/backlog/{priority}-{slug}.md`) or from the backlog item body. The field format is a bullet list following the header:

   ```markdown
   **Acceptance Criteria**:
   - {criterion 1}
   - {criterion 2}
   - {criterion 3}
   ```

   Parse each `-` line as a separate criterion.

5. Spawn a verification agent with subagent_type="general-purpose". Prompt must include: item title, plan file path, checklist status (100%), and each criterion listed individually as "Criterion N: {text}". Instruct the agent to: read the plan file, search `git log --oneline -20`, check relevant files for each criterion, and return per-criterion PASS/FAIL with file:line evidence. Required return format:

   <eg>
   [PASS] {criterion} — verified at {file}:{line} (or commit {sha})
   [FAIL] {criterion} — {gap description}
   Overall: PASS or FAIL (N/M criteria met)
   </eg>

   **If no acceptance criteria exist**: warn "No **Acceptance Criteria**: field found — falling back to description-based verification" and spawn agent with the description as the goal instead.

6. Parse the agent verdict:

   <eg>
   Acceptance Criteria Verification:

     [PASS] running X produces Y — verified at src/main.py:42
     [FAIL] file Z exists with field W — file exists but field W not found
     [PASS] test suite passes — confirmed via git log (commit abc123f)

   Overall: FAIL (2/3 criteria met)
   </eg>

7. Collect agent verdict:
   - **Overall PASS** (all criteria met): proceed to 9e
   - **Overall FAIL** (any criterion failed): report gaps, do not close:

<eg>
Verification FAILED for "{title}".

Per-criterion results:
{agent findings}

Address the failing criteria before closing.
</eg>

Then stop.

## 9e: Invoke backlog close

8. Invoke the backlog script (script updates per-item file and closes GitHub issue):

```bash
uv run .claude/skills/backlog/scripts/backlog.py close "{title}" --plan "{plan file path}" --checklist-pass -R Jamie-BitFlight/claude_skills
```

If invoked as `close #N`, use `#N` as the selector:

```bash
uv run .claude/skills/backlog/scripts/backlog.py close "#{N}" --plan "{plan file path}" --checklist-pass -R Jamie-BitFlight/claude_skills
```

9. Report the script output to the user.
