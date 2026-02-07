---
description: "List all git worktrees with enhanced status information"
allowed-tools: Bash, Read, TodoWrite
---

# Git Worktree List

Display all worktrees with detailed status information.

## Phase 1: Gather Worktree Information

**Get basic list**:
```bash
git worktree list --porcelain
```

**For each worktree, collect**:
1. Path and branch name
2. HEAD commit (short hash + message)
3. Worktree status (clean/dirty)

## Phase 2: Check Each Worktree Status

**For each worktree path**:

```bash
# Uncommitted changes count
git -C <path> status --porcelain | wc -l

# Unpushed commits (if remote exists)
git -C <path> rev-list --count @{upstream}..HEAD 2>/dev/null || echo "no-remote"

# Behind upstream
git -C <path> rev-list --count HEAD..@{upstream} 2>/dev/null || echo "0"
```

## Phase 3: Detect Stale Worktrees

**Check for stale entries**:
```bash
git worktree list --porcelain | grep "^worktree" | cut -d' ' -f2-
```

**For each path**:
- If directory doesn't exist → mark as **stale**
- If `.git` file is corrupted → mark as **stale**

**If stale worktrees found**:
- List them separately
- Suggest running `git worktree prune`

## Phase 4: Display Results

**Output format (table)**:

```
Git Worktrees:
┌─────────────────────────────────┬──────────────┬────────┬──────────┬──────────┐
│ Path                            │ Branch       │ Status │ Unpushed │ Behind   │
├─────────────────────────────────┼──────────────┼────────┼──────────┼──────────┤
│ /path/to/main                   │ main         │ clean  │ 0        │ 0        │
│ /path/to/feature-login          │ feature/login│ dirty  │ 3        │ 2        │
│ /path/to/bugfix                 │ bugfix/123   │ clean  │ 1        │ 0        │
└─────────────────────────────────┴──────────────┴────────┴──────────┴──────────┘
```

**Status indicators**:
| Status | Meaning |
|--------|---------|
| clean | No uncommitted changes |
| dirty | Has uncommitted changes |
| stale | Directory missing |

**Color hints** (in description):
- dirty + unpushed > 0 → needs attention
- behind > 0 → suggest pull/rebase

## Phase 5: Summary & Suggestions

**Show summary**:
- Total worktrees: N
- Clean: X, Dirty: Y, Stale: Z

**Suggestions based on status**:
| Condition | Suggestion |
|-----------|------------|
| Has stale worktrees | Run `git worktree prune` |
| Has dirty worktrees | Commit or stash changes |
| Has unpushed commits | Push changes |
| Branch merged to main | Consider removing worktree |

## Quick Actions

After listing, offer quick actions:
1. Prune stale worktrees
2. Show details for a specific worktree
3. Remove a worktree (`/git:worktree-remove`)

## Output Language

- All messages: Chinese
- Paths and branches: As-is
