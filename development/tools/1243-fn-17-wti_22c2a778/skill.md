# fn-17-wti Centralize Review Backend Detection

## Problem

Currently, 6 skill files detect `rp-cli` and `codex` availability at runtime using `which rp-cli` and `which codex`. This causes:

1. **LLM deviation** - Agents sometimes check wrong binary names (`rp`, `repoprompt` instead of `rp-cli`)
2. **Redundant checks** - Same detection repeated in every skill invocation (12+ subprocess calls)
3. **Inconsistency** - Ralph mode already handles this correctly via config

## Solution

Move detection to `/flow-next:setup` (one-time), persist to `.flow/config.json`, remove runtime detection from skills.

## Design Decisions

| Question | Decision |
|----------|----------|
| Tool unavailable at runtime but configured | Error with clear message: "rp-cli configured but not found. Run /flow-next:setup to reconfigure." |
| Setup when no tools available | Show options with note "neither detected", default to "none" |
| Config key | `review.backend` (nested, consistent with `memory.enabled`, `planSync.enabled`) |
| Re-run setup behavior | Pre-select current value, allow user to change |
| User skips question | Write `"none"` explicitly (no ambiguity) |

## Priority Order (all skills)

Skills use this priority (first match wins):
1. `--review=rp|codex|export|none` argument (explicit user choice)
2. `FLOW_REVIEW_BACKEND` env var (Ralph mode sets this)
3. `.flow/config.json` → `review.backend` (setup configured this)
4. **ERROR** - "No review backend configured. Run /flow-next:setup or pass --review=X"

**NO auto-detection fallback.** Trust explicit config or fail fast.

## Scope

### In scope
- Add review backend question to `/flow-next:setup`
- Add `review.backend` to default config in flowctl.py
- Remove runtime detection from 6 skill files
- Update error messages to guide users

### Out of scope
- Changes to `/flow-next:ralph-init` (already works correctly via `FLOW_REVIEW_BACKEND` env)
- Per-command backend config (future enhancement)
- Tool path caching (future enhancement)

## Files to Modify

| File | Change |
|------|--------|
| `plugins/flow-next/scripts/flowctl.py` | Add `review.backend` to `get_default_config()` |
| `plugins/flow-next/skills/flow-next-setup/workflow.md` | Add review backend question to Step 6 |
| `plugins/flow-next/skills/flow-next-plan/SKILL.md` | Remove `which` detection |
| `plugins/flow-next/skills/flow-next-work/SKILL.md` | Remove `which` detection |
| `plugins/flow-next/skills/flow-next-plan-review/SKILL.md` | Remove `which` detection |
| `plugins/flow-next/skills/flow-next-plan-review/workflow.md` | Remove `which` detection |
| `plugins/flow-next/skills/flow-next-impl-review/SKILL.md` | Remove `which` detection |
| `plugins/flow-next/skills/flow-next-impl-review/workflow.md` | Remove `which` detection |

## Quick Commands

```bash
# Run smoke tests
plugins/flow-next/scripts/smoke_test.sh

# Verify config default includes review.backend
python3 plugins/flow-next/scripts/flowctl.py init --json
cat .flow/config.json | jq '.review'

# Verify no 'which rp-cli' in skill files
grep -r "which rp-cli" plugins/flow-next/skills/ && echo "FAIL: detection still present" || echo "PASS: detection removed"
```

## Acceptance

- [ ] `/flow-next:setup` asks about review backend and writes to config
- [ ] `flowctl init` creates config with `review.backend: null` (or absent)
- [ ] All 6 skill files no longer run `which rp-cli` or `which codex`
- [ ] Skills read from: flag → env → config → error (no auto-detect)
- [ ] Error message is clear when no backend configured
- [ ] Ralph mode unaffected (FLOW_REVIEW_BACKEND env still works)
- [ ] Smoke tests pass
- [ ] README docs updated (remove "auto-detect" from priority)

## References

- Gap analysis identified edge cases and priority questions
- Ralph mode uses `FLOW_REVIEW_BACKEND` env (lines 798, 816 of ralph.sh)
- Config pattern in flowctl.py lines 76-78, 107-149
