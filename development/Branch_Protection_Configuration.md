---
name: Branch Protection Configuration
source: https://raw.githubusercontent.com/LearningCircuit/local-deep-research/main/.github/BRANCH_PROTECTION.md
original_path: .github/BRANCH_PROTECTION.md
source_repo: LearningCircuit/local-deep-research
category: development
subcategory: coding
tags: ['development']
collected_at: 2026-01-31T19:32:40.646780
file_hash: 6ee70b182c61aa2317d19d723b2f9e8dde0e28dfc99b899158fbba406a459684
---

# Branch Protection Configuration

Required settings for OSSF Scorecard compliance on `main` and `dev` branches.

## Required Settings (GitHub UI: Settings -> Branches -> Add rule)

### Protect matching branches
- [x] Require a pull request before merging
  - [x] Require at least 1 approving review
  - [x] Dismiss stale PR approvals when new commits are pushed
- [x] Require status checks to pass before merging
  - [x] Require branches to be up to date before merging
- [x] Do not allow bypassing the above settings

### Status checks to require
- `test` (unit tests)
- `lint` (code quality)
- `typecheck` (mypy)

## Setup Instructions

1. Go to repository **Settings** -> **Branches**
2. Click **Add branch protection rule**
3. Set "Branch name pattern" to `main`
4. Configure the settings listed above
5. Click **Create**
6. Repeat for `dev` branch

## Why This Matters

OSSF Scorecard's Branch-Protection check verifies that:
- Direct pushes to protected branches are prevented
- Code review is required before merging
- Status checks must pass before merging

These settings help prevent accidental or malicious changes to production code.

## References

- [OSSF Scorecard Branch-Protection Check](https://github.com/ossf/scorecard/blob/main/docs/checks.md#branch-protection)
- [GitHub Branch Protection Documentation](https://docs.github.com/en/repositories/configuring-branches-and-merges-in-your-repository/managing-protected-branches/about-protected-branches)
