---
description: "Remove a git worktree and optionally delete the branch"
argument-hint: Worktree path or branch name, or "--merged" for batch cleanup, or "--prune" for stale cleanup
allowed-tools: Bash, Read, TodoWrite, AskUserQuestion
---

# Git Worktree Remove

Safely remove a worktree and clean up associated resources.

## Mode Detection

**Parse `$ARGUMENTS`**:
| Argument | Mode |
|----------|------|
| `--merged` | Batch remove all merged worktrees |
| `--prune` | Auto-detect and clean stale worktrees |
| `--all` | Interactive batch removal |
| path/branch | Single worktree removal |
| (empty) | Interactive selection |

---

## Mode: Stale Cleanup (`--prune`)

### Step 1: Detect Stale Worktrees

```bash
# Get all registered worktree paths
git worktree list --porcelain | grep "^worktree" | cut -d' ' -f2-
```

**For each path, check**:
```bash
[ -d "<path>" ] && echo "exists" || echo "stale"
```

### Step 2: Show Stale Worktrees

**If stale found**:
```
发现 N 个 stale worktree（目录已不存在）：
- /path/to/deleted-worktree-1
- /path/to/deleted-worktree-2
```

**If none found**:
```
没有发现 stale worktree，所有 worktree 目录都存在。
```

### Step 3: Prune

**Ask user to confirm**, then:
```bash
git worktree prune -v
```

**Show result**: List of pruned entries.

---

## Mode: Batch Merged Cleanup (`--merged`)

### Step 1: Find Merged Worktrees

```bash
# Get merged branches
git branch --merged main | grep -v "main\|master\|\*"
```

**Cross-reference with worktrees**:
```bash
git worktree list --porcelain
```

**Build list**: Worktrees whose branches are already merged to main.

### Step 2: Show Candidates

```
发现 N 个已合并的 worktree 可以清理：
┌─────────────────────────────────┬──────────────────┬────────────────┐
│ Path                            │ Branch           │ Merged To      │
├─────────────────────────────────┼──────────────────┼────────────────┤
│ ../project-feature-a            │ feature/a        │ main           │
│ ../project-bugfix-123           │ bugfix/123       │ main           │
└─────────────────────────────────┴──────────────────┴────────────────┘
```

### Step 3: Confirm Batch Removal

**Ask user**:
- Remove all? (Y/n)
- Also delete branches? (Y/n)
- Also delete remote branches? (y/N)

### Step 4: Execute Batch Removal

**For each worktree**:
```bash
git worktree remove <path>
git branch -d <branch>  # if user confirmed
git push origin --delete <branch>  # if user confirmed
```

**Show progress**: `[1/N] 移除 feature/a ...`

### Step 5: Final Cleanup

```bash
git worktree prune
git fetch --prune  # if remote branches deleted
```

---

## Mode: Interactive Batch (`--all`)

### Step 1: List All Worktrees with Status

Show table with checkboxes:
```
选择要删除的 worktree：
[ ] ../project-feature-a     feature/a     dirty   (3 uncommitted)
[x] ../project-bugfix-123    bugfix/123    clean   (merged)
[ ] ../project-experiment    experiment    clean   (not merged)
```

### Step 2: User Selection

Ask user to select worktrees to remove (comma-separated indices or "merged" for all merged).

### Step 3: Safety Checks

**For each selected**:
- If dirty: warn and require explicit confirm
- If not merged: warn and require explicit confirm

### Step 4: Execute

Same as batch removal with progress indicator.

---

## Mode: Single Removal (Default)

## Phase 1: List Worktrees

**Actions**:
```bash
git worktree list
```

**Show table**:
| Path | Branch | Status |
|------|--------|--------|
| /path/to/main | main | (bare) |
| /path/to/feature | feature/x | clean |

## Phase 2: Select Worktree

**If `$ARGUMENTS` provided**:
- Match by path or branch name
- Validate worktree exists

**If empty**:
- Show worktree list
- Ask user to select one

**Prevent removing main worktree**.

## Phase 3: Pre-Removal Checks

**Check for uncommitted changes**:
```bash
cd <worktree-path> && git status
```

**If changes exist**:
- Show changed files
- Ask user:
  - Stash changes?
  - Commit changes?
  - Discard changes?
  - Abort removal?

**Check for unpushed commits**:
```bash
git log origin/<branch>..<branch> --oneline
```

If unpushed: warn user and ask to push first.

## Phase 4: Remove Worktree

**Standard removal**:
```bash
git worktree remove <path>
```

**Force removal** (if dirty):
```bash
git worktree remove --force <path>
```

## Phase 5: Branch Cleanup

**Ask user**: Delete the branch too?

**If yes, check merge status**:
```bash
git branch --merged main | grep <branch>
```

**Delete branch**:
- If merged: `git branch -d <branch>`
- If not merged: warn, require `git branch -D <branch>`

**Delete remote branch** (optional):
```bash
git push origin --delete <branch>
```

## Phase 6: Prune Stale References

**Clean up**:
```bash
git worktree prune
```

**Show final worktree list**.

## Safety Checks

| Check | Action |
|-------|--------|
| Uncommitted changes | Require user decision |
| Unpushed commits | Warn and confirm |
| Unmerged branch | Require force confirm |
| Main worktree | Block removal |

## Output Language

- All messages: Chinese
- Paths and commands: As-is
