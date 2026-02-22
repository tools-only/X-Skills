---
name: clean-branches
description: >
  This skill should be used when the user says "clean up branches", "delete merged
  branches", "prune stale branches", "git branch cleanup", "remove old branches",
  or wants to tidy up or purge old branches.
model: haiku
argument-hint: [branch-pattern] - optional pattern to filter branches
allowed-tools:
  - Bash
  - AskUserQuestion
---

# Clean Git Branches

Safely remove merged and stale git branches with confirmation.

## Process

**0. Parse arguments**
If `$ARGUMENTS` provided, treat it as a glob pattern to filter branch candidates to only matching branches (e.g., `feature/*` shows only feature branches). Applied in Step 2 to both merged and stale detection.

**1. Fetch latest state**
```bash
git fetch --all --prune
```
If git fetch fails (no remotes configured), note that remote branch data is unavailable and continue with local analysis only.

**2. Identify candidates**

Find merged branches (apply pattern filter if provided):
```bash
MERGED=$(git branch --merged main | grep -v "^\*\|main\|master\|develop")
[ -n "$ARGUMENTS" ] && MERGED=$(echo "$MERGED" | grep "$ARGUMENTS")
echo "$MERGED"
```

Find stale branches (no commits in 30+ days, apply pattern filter if provided):
```bash
CUTOFF=$(python3 -c "import time; print(int(time.time()) - 30*86400)")
git for-each-ref --sort=-committerdate \
  --format='%(refname:short) %(committerdate:unix) %(committerdate:relative)' \
  refs/heads/ | while read branch ts reldate; do
  [ -n "$ARGUMENTS" ] && [[ "$branch" != $ARGUMENTS ]] && continue
  if (( ts < CUTOFF )); then
    echo "$branch ($reldate)"
  fi
done
```
Unix timestamps are used for accurate threshold comparison — git's relative date strings ("5 weeks ago") would miss branches 35–59 days old if matched by pattern.

**3. Categorize and display**

Present branches in three groups: merged (safe to delete), stale (no recent commits), and protected (never touch). Show branch name and age for stale entries.

```
MERGED BRANCHES (safe to delete):
- feature/old-feature
- fix/completed-fix

STALE BRANCHES (no recent commits):
- experiment/abandoned (3 months ago)
- wip/forgotten (6 months ago)

PROTECTED (never delete):
- main, master, develop
```

If both lists are empty, report "No branches to clean" and stop.

**4. Confirm before deletion**

Use AskUserQuestion:
- "Delete all merged branches?"
- "Delete specific stale branches?" (list options)
- "Skip and keep all?"

**5. Execute deletion**

Local only (safe):
```bash
git branch -d <branch-name>
```

If user explicitly requests remote cleanup:
```bash
git push origin --delete <branch-name>
```

## Safety Rules

- Never delete: main, master, develop, release/*
- Always confirm before any deletion
- Use `-d` not `-D`: `-d` refuses to delete branches with unmerged commits, so git itself enforces the safety check; `-D` bypasses it
- Show what will be deleted before acting
- Remote deletion requires explicit user confirmation; never delete remotes unless the user says so directly

## Output

Summary:
- Branches deleted (local)
- Branches deleted (remote, if requested)
- Branches kept
- Any errors encountered
