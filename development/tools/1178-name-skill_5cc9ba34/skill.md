---
name: push-pr
description: >
  This skill should be used when the user wants to push commits, create or update a pull
  request, or submit code for review. Triggers on "push this", "push my changes",
  "create a PR", "open a pull request", "make a PR", "submit for review", "send this up",
  "open PR", "pr please", or any mention of pushing code or creating a PR.
model: sonnet
context: fork
argument-hint: "[status: 1=open|2=draft|3=ready] [base-branch]"
allowed-tools:
  - Bash(git:*)
  - Bash(gh:*)
  - Skill
  - AskUserQuestion
---

# Push & PR

Push commits and create/update pull requests with automatic branch management and
scope-aware multi-PR splitting.

## Arguments

Parse flexibly from `$ARGUMENTS`:
- **status**: `1`=opened, `2`=draft, `3`=ready (default: new PR=opened, update=draft)
- **base-branch**: Target branch (default: `main`)

## Workflow

### 1. Pre-Flight

```bash
git status --porcelain
git fetch origin
```

If uncommitted changes are detected, invoke `Skill: commit` to commit first.

### 2. Branch Management

```bash
CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
UNPUSHED=$(git rev-list origin/$CURRENT_BRANCH..HEAD --count 2>/dev/null || echo "0")
```

**If on main/master with unpushed commits — cut a feature branch before proceeding.**

Branch naming rule: prefix from the primary commit's conventional type (`feat/`, `fix/`,
`docs/`, `chore/`, `refactor/`); slug from the commit scope or subject, lowercase hyphens
only, max 45 chars total. Use the scope if present (`feat(auth)` → `feat/auth`); otherwise
condense the subject to 2–4 words (`fix login redirect timeout` → `fix/login-redirect`).

```bash
git log origin/main..HEAD --oneline
git checkout -b <derived-branch-name>
git checkout main && git reset --hard origin/main
git checkout <derived-branch-name>
```

**If already on a feature branch:** skip, proceed to step 3.

### 3. Context Gathering

```bash
BASE=${BASE_BRANCH:-main}
git log $BASE..HEAD --oneline --reverse
git diff $BASE...HEAD --stat
```

Record: commit count, conventional-commit types and scopes present, total diff lines
(approximate from `--stat` output).

### 4. Scope Analysis

Evaluate whether the changeset warrants multiple PRs. A split is warranted when either:
- **Size**: total diff exceeds ~400 lines (code lines; ignore lock files and generated files)
- **Diversity**: commits span 3+ distinct conventional-commit scopes or types (e.g.,
  `feat(auth)`, `fix(ui)`, `chore(deps)` are three distinct scopes)

**If neither condition is met:** proceed to step 5 as a single PR.

**If either condition is met — propose stacked PRs:**

Cluster commits by scope/type in the order they were made. Each cluster becomes one PR
targeting the previous cluster's branch (the first targets `$BASE`). Present the plan:

```
Proposed stacked PRs (each PR targets the previous branch):
  PR 1  [base: main]         feat/auth          — commits: abc1234, def5678
  PR 2  [base: feat/auth]    fix/ui-redirect    — commits: ghi9012
  PR 3  [base: fix/ui-...]   chore/cleanup      — commits: jkl3456
```

Use AskUserQuestion: "Split into N stacked PRs as shown above, or push as a single PR?"

**If user declines split:** proceed to step 5 as a single PR.

**If user confirms split — stacked PR execution:**

For each cluster in order:
1. Create a branch from the previous cluster's branch (`$BASE` for cluster 1):
   ```bash
   git checkout -b <cluster-branch> <previous-branch>
   git cherry-pick <sha1> <sha2> ...
   ```
2. Push: `git push -u origin <cluster-branch>`
3. Generate PR body (step 7) with `--base <previous-branch>`.
4. Create PR targeting the correct base branch.
5. Repeat for next cluster.

After all PRs are created, check out the last cluster's branch and report the full
stack (see Output). Exit — skip steps 5–7.

### 5. PR Status

```bash
BRANCH=$(git rev-parse --abbrev-ref HEAD)
gh pr list --head "$BRANCH" --json number,state
```

Use the provided status argument, or default: new PR=opened, update=draft.

### 6. Push

```bash
BRANCH=$(git rev-parse --abbrev-ref HEAD)
if git rev-parse --abbrev-ref --symbolic-full-name @{u} 2>/dev/null; then
    git push
else
    git push -u origin "$BRANCH"
fi
```

### 7. PR Creation/Update

Generate the PR body using the format script. Pass the correct base branch — for stacked
PRs this is the previous cluster's branch, not necessarily `$BASE`:

```bash
python3 ${CLAUDE_PLUGIN_ROOT}/skills/push-pr/scripts/format-pr-body.py --base "$BASE"
```

Exit 1 means no changes found relative to base; report to user. On success, use stdout
as the PR body directly.

**New PR:**
```bash
gh pr create --title "<title>" --body "<format-pr-body output>" --base "$BASE"
# If status=ready: gh pr ready
```

**Existing PR:** Add a comment listing new commits since last push; update PR status if
the status argument changed.

```bash
PR_NUM=$(gh pr list --head "$BRANCH" --json number -q '.[0].number')
gh pr comment $PR_NUM --body "New commits: ..."
```

## Constraints

- NO Co-authored-by or AI signatures
- NO "Generated with Claude Code"
- NO emojis in PR title or description
- Use existing git user config only

## Edge Cases

- No remote → suggest `git remote add origin <url>` and stop
- No `gh` CLI → report requirement and stop
- Branch behind remote → pull/rebase before pushing
- No commits to push → report and stop
- Cherry-pick conflict during stacked flow → stop, report which cluster failed, suggest
  resolving manually then re-running

## Output

Single PR: branch name, PR URL, PR status (opened/draft/ready).

Stacked PRs: ordered list showing each PR URL and the branch it targets, plus the name
of the final branch now checked out.

## Scripts

`scripts/format-pr-body.py` — generates a formatted markdown PR body by reading git
history and diff stats. Invoked in step 7.

```bash
# Usage
python3 ${CLAUDE_PLUGIN_ROOT}/skills/push-pr/scripts/format-pr-body.py --base main
python3 ${CLAUDE_PLUGIN_ROOT}/skills/push-pr/scripts/format-pr-body.py --base feat/auth --output json
```

Exit 0: PR body written to stdout. Exit 1: no changes found (script reports to stderr).
JSON mode returns `{"title": "<first commit message>", "body": "<markdown body>"}`.
