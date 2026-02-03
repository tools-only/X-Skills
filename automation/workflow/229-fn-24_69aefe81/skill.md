# fn-2.4 Update plan-review skill for codex backend

## Description

Update `/flow-next:plan-review` skill to support codex backend alongside rp.

### Changes needed

Same pattern as fn-2.3 (impl-review):

1. **workflow.md** - Add backend detection and branching
2. **Backend dispatch**:
   ```
   If backend == "codex":
       eval "$(flowctl codex plan-review $EPIC_ID)"
   ```

### Files to modify

- `plugins/flow-next/skills/flow-next-plan-review/workflow.md`
- `plugins/flow-next/skills/flow-next-plan-review/SKILL.md` (if needed)
## Acceptance
- [ ] Skill detects backend from env var or config
- [ ] Codex path calls `flowctl codex plan-review`
- [ ] RP path unchanged (no regression)
- [ ] Verdict returned correctly from both backends
## Done summary
- Updated SKILL.md with backend selection logic
- Updated workflow.md with codex backend workflow (Phase 0 + codex section)
- Labeled RP phases for clarity
- Updated anti-patterns for both backends

Why:
- Enable codex as alternative to RepoPrompt for plan reviews
- No regression on RP workflow

Verification:
- smoke_test.sh passes
## Evidence
- Commits: 4d2b1a5945ec615f37b6b0dbf1b59cfebe729f95
- Tests: smoke_test.sh
- PRs: