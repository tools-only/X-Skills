# fn-2.5 Update Ralph templates for codex option

## Description

Update Ralph templates to support codex as review backend option.

### Files to update

1. **config.env** - Add codex to options:
   ```bash
   # Review backends: rp (RepoPrompt, macOS), codex (cross-platform), none
   PLAN_REVIEW={{PLAN_REVIEW}}    # rp|codex|none
   WORK_REVIEW={{WORK_REVIEW}}    # rp|codex|none
   ```

2. **prompt_plan.md** - Add codex dispatch:
   ```markdown
   If PLAN_REVIEW=codex:
       eval "$(flowctl codex plan-review $EPIC_ID)"
   ```

3. **prompt_work.md** - Add codex dispatch:
   ```markdown
   If WORK_REVIEW=codex:
       eval "$(flowctl codex impl-review $TASK_ID --base $BASE_BRANCH)"
   ```

4. **ralph.sh** - Update UI to show "Codex" instead of "RepoPrompt" when applicable

### Files to modify

- `plugins/flow-next/skills/flow-next-ralph-init/templates/config.env`
- `plugins/flow-next/skills/flow-next-ralph-init/templates/prompt_plan.md`
- `plugins/flow-next/skills/flow-next-ralph-init/templates/prompt_work.md`
- `plugins/flow-next/skills/flow-next-ralph-init/templates/ralph.sh`
## Acceptance
- [ ] config.env documents `rp|codex|none` options
- [ ] prompt_plan.md dispatches to codex when configured
- [ ] prompt_work.md dispatches to codex when configured
- [ ] ralph.sh UI shows "Codex" when PLAN_REVIEW=codex or WORK_REVIEW=codex
- [ ] Ralph loop works with PLAN_REVIEW=codex WORK_REVIEW=codex
## Done summary
Updated Ralph templates for codex backend support: config.env with option docs, prompt_plan.md and prompt_work.md with codex dispatch, ralph.sh with UI and receipt verification for both rp and codex modes.
## Evidence
- Commits: 72263cd
- Tests: smoke_test.sh
- PRs: