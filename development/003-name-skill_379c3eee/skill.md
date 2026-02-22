---
name: commit
description: >
  This skill should be used when the user says "commit my changes", "commit this",
  "git commit", "save my work", "stage and commit", or mentions committing code.
  Also triggers on "create a commit" or "commit what I've done".
model: sonnet
context: fork
agent: general-purpose
argument-hint: "[push]"
allowed-tools:
  - Bash(git:*)
  - Bash(cargo:*)
  - Bash(npm:*)
  - Bash(ruff:*)
  - AskUserQuestion
---

# Commit

Analyze uncommitted changes and create well-organized commits using conventional commit format.

## Workflow

### 1. Discover Changes

```bash
git status --porcelain
git diff --stat
```

If no changes, report "Nothing to commit" and stop.

### 2. Stage Files

Run `git add -A` to stage all changes.

**Exclude temporary files** matching: `scratch.*`, `temp.*`, `debug.*`, `playground.*`, `*.log`, `dist/`, `build/`, `target/`, `node_modules/`.

If temporary files detected:
1. Unstage with `git reset HEAD <file>`
2. Ask user if they want to add to .gitignore

### 3. Analyze Commit Boundaries

For each changed file, write a one-line PURPOSE description (not file location).

**Group by PURPOSE, not directory:**
- Same goal = one commit
- Different goals = separate commits

Each commit should represent one logical change because atomic commits enable `git bisect` and `git revert` without side-effects.

**Signs of separate concerns:**
- "Added X" AND "Fixed Y" (feature + bugfix)
- Changes that could be reverted independently

If multiple concerns: use `git reset HEAD` then `git add <specific-files>` for each group. Commit foundational changes first.

**Handle renames (R status):** When splitting, add BOTH old and new paths. Git detects renames by similarity scoring across the old/new pair — staging only the new path causes git to log a delete + add, losing rename history.

### 4. Validate (if available)

Run the validation script after staging:

```bash
python3 ${CLAUDE_PLUGIN_ROOT}/skills/commit/scripts/validate.py . --output json
```

Interpret the result:
- Exit 0 → validation passed; proceed to Step 5
- Exit 1 → validation failed; parse the `output` field, report the error to the user, and stop
- Exit 2 → no validator found for this project; skip gracefully and proceed to Step 5

### 5. Create Commit

Check recent commit style:
```bash
git log --oneline -10
```

Use conventional commit format:
```
<type>(<scope>): <description>
```

**Types:** feat, fix, docs, refactor, test, chore, perf

**Rules:**
- Lowercase, no period, imperative mood
- Max 72 chars for subject
- NO Co-authored-by trailers
- NO AI attribution
- NO emojis

### 6. Push (if requested)

If user mentions "push" or arguments contain "push", run `git push`. If push fails, report the error and stop — do not retry or force-push.

## Output

One line per commit (hash + message). If temporary files were excluded, list them as bullets below.

## Scripts

`scripts/validate.py` — detects project type and runs the appropriate linter. Invoked in Step 4.

```bash
# Usage
python3 ${CLAUDE_PLUGIN_ROOT}/skills/commit/scripts/validate.py .
python3 ${CLAUDE_PLUGIN_ROOT}/skills/commit/scripts/validate.py /path/to/project --output json
```
