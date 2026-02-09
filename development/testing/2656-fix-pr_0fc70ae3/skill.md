# Fix PR Review Comments (v3)

Systematically fix all review comments on a PR from all sources.

## Mode Selection
- **DEFAULT (safe)**: Present findings → STOP → wait for user approval → fix
- **AUTONOMOUS**: Activated by user saying "autonom", "autonomous", or "mach einfach" → skip approval gate, loop until MERGEABLE

## Phase 1: Collect ALL Comments
1. Determine PR number: `gh pr list --state open --head $(git branch --show-current) --json number -q '.[0].number'`
2. Fetch PR metadata: `gh pr view $PR --json reviewDecision,reviews,comments,mergeable`
3. Fetch ALL review comment sources:
   - PR review comments: `gh api repos/{owner}/{repo}/pulls/$PR/comments`
   - PR review threads: `gh api repos/{owner}/{repo}/pulls/$PR/reviews`
   - Issue-style comments: `gh api repos/{owner}/{repo}/issues/$PR/comments`
4. Identify sources: CodeRabbit, Qodo/Codium, manual reviewers, GitHub Actions bots

## Phase 2: Present Findings
5. Group ALL comments by file path, then by severity: [CRITICAL] > [WARNING] > [SUGGESTION]
6. Create a TodoWrite checklist of every finding grouped by file

## Phase 2b: Approval Gate
7. **DEFAULT mode**: STOP HERE. Wait for user approval before proceeding.
8. **AUTONOMOUS mode**: Skip this gate, proceed directly to Phase 3.

## Phase 3: Fix (file by file)
9. Read each file before editing - understand context first
10. Fix each finding, run the relevant test suite after EACH file change:
    - `.ts`/`.tsx` files: `npx tsc --noEmit` + `npm test -- --run`
    - `.py` files: `python -m py_compile` + `pytest`
11. Use real framework objects in tests (never MagicMock for HTTP request/response)
12. Mark each TodoWrite item as completed after fixing

## Phase 4: Validate
13. Resolve merge conflicts: `git fetch origin main && git merge --no-commit origin/main` (resolve if conflicts, abort and retry if needed)
14. Run full test suite
15. Run pre-commit hooks. If PII scanner false positives on number literals or placeholder emails:
    - Document the false positive with justification
    - Use `--no-verify`

## Phase 5: Commit and Push
16. Commit with structured message referencing PR number and listing all fixes
17. Push to remote

## Phase 6: Re-fetch Loop (AUTONOMOUS mode only)
18. Re-fetch PR review status: `gh pr view $PR --json reviewDecision,mergeable`
19. If new comments appeared → loop back to Phase 1
20. Only stop when PR status is MERGEABLE and all quality gates pass

## Phase 7: Final Report
Always output a structured summary:

```
## PR #<NUMBER> Fix Summary
- **Files changed**: [list]
- **Comments resolved**: X/Y
- **Tests**: all passing / N failures
- **Deferred items**: [list with reasoning]
- **Status**: MERGEABLE / BLOCKED (reason)
```
