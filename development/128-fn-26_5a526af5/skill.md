# fn-2.6 Update ralph-guard.py for codex patterns

## Description

Update ralph-guard.py hooks to validate codex calls same as rp calls.

### Changes needed

1. **Block direct codex calls** (prevent drift):
   ```python
   BLOCKED_CODEX_PATTERNS = [
       r'codex\s+exec\s+(?!.*--sandbox)',      # exec without sandbox
       r'codex\s+review',                       # direct review (should use flowctl)
       r'codex\s+exec\s+.*--last',             # --last flag (use explicit session_id)
   ]
   ```

2. **Allow flowctl wrappers**:
   ```python
   ALLOWED_PATTERNS = [
       r'flowctl\s+codex\s+',    # All flowctl codex commands OK
       r'flowctl\s+rp\s+',       # All flowctl rp commands OK
   ]
   ```

3. **Track codex mode** - Similar state tracking as rp mode

4. **Session enforcement** (optional, warn-only):
   - If receipt exists with `session_id` and codex exec starts new session → warn
   - Block `--last` flag entirely (must use explicit session_id)

### Files to modify

- `plugins/flow-next/scripts/hooks/ralph-guard.py`
## Acceptance
- [ ] Direct `codex exec` without sandbox blocked
- [ ] Direct `codex review` blocked (must use flowctl wrapper)
- [ ] `codex exec --last` blocked (must use explicit session_id)
- [ ] `flowctl codex *` commands allowed
- [ ] Existing `flowctl rp *` validation unchanged
- [ ] State tracking works for codex mode
## Done summary
Updated ralph-guard.py to support codex backend: tracks codex review state, validates receipt writes for both backends, detects verdict from codex output.
## Evidence
- Commits: a0c302a
- Tests: smoke_test.sh
- PRs: