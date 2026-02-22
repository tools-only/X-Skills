---
skill_name: git-workflows
skill_category: devops
description: Git workflows, branching strategies, commit conventions, and version control best practices
allowed_tools: [Bash, Read, Edit, Glob, Grep]
token_estimate: 1200
version: 1.0
last_updated: 2026-01-21
owner: DevOps Team
status: active

tags: [git, version-control, branching, commits, workflow, collaboration]
related_skills: [ci-cd-patterns, github-actions]

trigger_files: [".git/config", ".gitignore", ".gitattributes", "CONTRIBUTING.md"]
trigger_keywords: [git, commit, branch, merge, rebase, pull request, PR, version control, gitflow, trunk-based]

quality_keywords: [conventional-commits, branching-strategy, merge-conflict, code-review]
---

# Git Workflows

Best practices for Git version control, branching strategies, and commit conventions.

## Purpose

- Establish consistent Git workflows across projects
- Ensure clean, meaningful commit history
- Prevent common Git mistakes and conflicts
- Facilitate effective code collaboration

---

## Core Operations

### Checking Status

```bash
# Full status
git status

# Short status
git status -s

# Show branch tracking info
git status -sb
```

### Staging and Committing

```bash
# Stage specific files
git add path/to/file.ts

# Stage all changes
git add .

# Stage interactively (review each change)
git add -p

# Commit with message
git commit -m "feat(scope): description"

# Commit with multi-line message
git commit -m "feat(scope): short description" \
           -m "- Detail 1" \
           -m "- Detail 2"
```

### Viewing History

```bash
# Recent commits (oneline)
git log --oneline -10

# Commits with graph
git log --oneline --graph --all

# Commits by author
git log --author="name" --oneline

# Diff of last commit
git diff HEAD~1
```

### Branching

```bash
# Create and switch to new branch
git checkout -b feature/my-feature

# Switch to existing branch
git checkout main

# List branches (local and remote)
git branch -a

# Delete merged branch
git branch -d feature/done

# Delete unmerged branch (force)
git branch -D feature/abandoned
```

### Syncing with Remote

```bash
# Fetch latest (no merge)
git fetch origin

# Pull with rebase (cleaner history)
git pull --rebase origin main

# Push to remote
git push origin feature/my-feature

# Push and set upstream
git push -u origin feature/my-feature
```

---

## Commit Message Convention

Use **Conventional Commits** format:

```
<type>(<scope>): <description>

[optional body]

[optional footer]
```

### Types

| Type | When to Use |
|------|-------------|
| `feat` | New feature for users |
| `fix` | Bug fix for users |
| `docs` | Documentation only |
| `style` | Formatting, no code change |
| `refactor` | Code restructuring, no behavior change |
| `perf` | Performance improvement |
| `test` | Adding/fixing tests |
| `chore` | Build, tooling, dependencies |

### Examples

```bash
# Feature
git commit -m "feat(auth): add OAuth2 login support"

# Bug fix
git commit -m "fix(api): handle null response from payment service"

# Refactor
git commit -m "refactor(db): extract connection pool to separate module"

# Breaking change
git commit -m "feat(api)!: change response format to JSON:API spec" \
           -m "BREAKING CHANGE: All endpoints now return JSON:API format"
```

---

## Branching Strategies

### Trunk-Based Development (Recommended)

**Best for:** Teams with strong CI/CD, frequent releases.

```
main ────●────●────●────●────●────●────
          \        /
           ●──────●  (short-lived feature branch)
```

**Rules:**
- Feature branches live < 2 days
- Merge to main frequently
- Use feature flags for incomplete work
- main is always deployable

### GitFlow

**Best for:** Scheduled releases, multiple versions in production.

```
main     ────●─────────────●─────────────●────
              \           / \           /
hotfix         ●─────────●   \         /
                              \       /
develop  ────●────●────●────●──●────●────●────
              \        /
feature        ●──────●
```

**Branches:**
- `main` - Production releases
- `develop` - Integration branch
- `feature/*` - New features
- `release/*` - Release preparation
- `hotfix/*` - Production fixes

---

## Common Workflows

### Starting New Feature

```bash
# Update main
git checkout main
git pull --rebase origin main

# Create feature branch
git checkout -b feature/user-authentication

# Work, commit, repeat
git add .
git commit -m "feat(auth): add login form component"

# Push to remote
git push -u origin feature/user-authentication
```

### Updating Feature Branch with Main

```bash
# Fetch latest
git fetch origin

# Rebase on main (cleaner than merge)
git rebase origin/main

# If conflicts, resolve then:
git add .
git rebase --continue

# Force push (required after rebase)
git push --force-with-lease
```

### Creating Pull Request

```bash
# Ensure branch is up to date
git fetch origin
git rebase origin/main

# Push final changes
git push origin feature/my-feature

# Create PR via CLI
gh pr create --title "feat(auth): add OAuth2 support" \
             --body "## Summary\n- Added OAuth2 flow\n- Added token refresh"
```

### Handling Merge Conflicts

```bash
# During rebase, if conflict occurs:
# 1. Edit conflicted files (look for <<<<<<< markers)
# 2. Stage resolved files
git add path/to/resolved-file.ts

# 3. Continue rebase
git rebase --continue

# Or abort if needed
git rebase --abort
```

---

## Anti-Patterns

### Anti-Pattern 1: Force Push to Shared Branches

| Aspect | Description |
|--------|-------------|
| **WHY** | Overwrites others' work. Causes sync issues. Can lose commits. |
| **DETECTION** | `git push --force` to main, develop, or any shared branch. |
| **FIX** | Use `--force-with-lease` on feature branches only. Never force push shared branches. |

### Anti-Pattern 2: Committing Secrets

| Aspect | Description |
|--------|-------------|
| **WHY** | Secrets in git history are permanent. Even after removal, they're in history. |
| **DETECTION** | API keys, passwords, tokens in committed files. |
| **FIX** | Use `.gitignore` for `.env` files. Use secret managers. If leaked, rotate credentials immediately. |

```bash
# Add to .gitignore BEFORE committing
echo ".env" >> .gitignore
echo "*.pem" >> .gitignore
echo "credentials.json" >> .gitignore
```

### Anti-Pattern 3: Giant Commits

| Aspect | Description |
|--------|-------------|
| **WHY** | Hard to review. Difficult to revert. Obscures history. |
| **DETECTION** | Commits touching 20+ files. Commits with "various changes" messages. |
| **FIX** | Commit early, commit often. One logical change per commit. |

### Anti-Pattern 4: Merging Without Updating

| Aspect | Description |
|--------|-------------|
| **WHY** | Creates unnecessary merge commits. Can introduce conflicts at merge time. |
| **DETECTION** | Feature branches far behind main. Merge conflicts discovered at PR time. |
| **FIX** | Rebase on main daily. Keep feature branches short-lived. |

---

## Quick Reference

| Task | Command |
|------|---------|
| Check status | `git status -sb` |
| Stage all | `git add .` |
| Commit | `git commit -m "type(scope): msg"` |
| New branch | `git checkout -b feature/name` |
| Update from main | `git rebase origin/main` |
| Push (first time) | `git push -u origin branch` |
| Push (subsequent) | `git push` |
| View log | `git log --oneline -10` |
| Undo last commit | `git reset --soft HEAD~1` |
| Discard changes | `git checkout -- file.ts` |

---

## Safety Rules

1. **Never force push to main/develop** - Use `--force-with-lease` only on feature branches
2. **Never commit secrets** - Use `.gitignore` and secret managers
3. **Never rewrite shared history** - Only rebase unpushed or personal branches
4. **Always pull before push** - Avoid unnecessary merge conflicts
5. **Keep commits atomic** - One logical change per commit
