# fn-2.8 Update docs for codex backend

## Description

Update all documentation to include Codex as cross-platform alternative to RepoPrompt.

**Key message**: RepoPrompt remains the recommended option (macOS), Codex is the cross-platform fallback for Windows/Linux users.

### Files to update

1. **plugins/flow-next/README.md**:
   - Add Codex to review backend options
   - Keep RepoPrompt as primary recommendation
   - Note Codex for Windows/Linux users

2. **plans/ralph-e2e-notes.md**:
   - Add Codex backend section
   - Document session continuity behavior
   - Add codex-specific troubleshooting

3. **plans/ralph-getting-started.md**:
   - Add Codex setup option alongside RP
   - Document `PLAN_REVIEW=codex WORK_REVIEW=codex` config

4. **CLAUDE.md**:
   - Update review backend options

### Tone

- "RepoPrompt (macOS, recommended)" - keep referral link
- "Codex (cross-platform alternative)" - for Windows/Linux users
- Don't oversell Codex - it's a fallback, RP is better UX on macOS

## Acceptance
- [ ] flow-next README documents both backends
- [ ] ralph-e2e-notes.md includes Codex section
- [ ] ralph-getting-started.md shows Codex setup
- [ ] RepoPrompt referral link preserved and recommended
- [ ] Codex positioned as cross-platform alternative, not replacement

## Done summary
Added codex backend documentation: flowctl.md with codex command group, ralph.md with integration section and config table, README.md with updated requirements and review sections.
## Evidence
- Commits: 1349589
- Tests: smoke_test.sh
- PRs: