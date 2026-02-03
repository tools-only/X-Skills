# fn-17-wti.7 Run smoke tests and update docs

## Description
Run smoke tests to verify changes, update README docs to remove "auto-detect" from priority.

## Steps

### 1. Run smoke tests

```bash
plugins/flow-next/scripts/smoke_test.sh
```

All tests must pass.

### 2. Verify no detection code remains

```bash
# Should return no matches
grep -r "which rp-cli" plugins/flow-next/skills/ && echo "FAIL" || echo "PASS: no which rp-cli"
grep -r "which codex" plugins/flow-next/skills/ && echo "FAIL" || echo "PASS: no which codex"
grep -r "HAVE_RP" plugins/flow-next/skills/ && echo "FAIL" || echo "PASS: no HAVE_RP"
grep -r "HAVE_CODEX" plugins/flow-next/skills/ && echo "FAIL" || echo "PASS: no HAVE_CODEX"
```

### 3. Test config default

```bash
rm -rf /tmp/test-config && mkdir /tmp/test-config && cd /tmp/test-config
python3 /path/to/flowctl.py init --json
cat .flow/config.json | jq '.review'
# Should output: {"backend": null}
```

### 4. Update README.md docs

In `plugins/flow-next/README.md`, find the Backend Selection / Review Backend section.

Update priority from:
```
1. `--review=X` flag
2. `FLOW_REVIEW_BACKEND` env
3. `.flow/config.json` → `review.backend`
4. Auto-detect (rp-cli preferred)
```

To:
```
1. `--review=X` flag
2. `FLOW_REVIEW_BACKEND` env
3. `.flow/config.json` → `review.backend`
4. Error: "No review backend configured"
```

Add note: "Run /flow-next:setup to configure your preferred review backend."

### 5. Update CHANGELOG.md

Add entry for version bump (if not already done):
```
- Remove runtime `which rp-cli`/`which codex` detection from skills
- Add review backend question to /flow-next:setup
- Skills now require explicit config (flag/env/setup) - no auto-detect
```

### 6. Version bump

Bump version in:
- `plugins/flow-next/.claude-plugin/plugin.json`
- `.claude-plugin/marketplace.json`
- `plugins/flow-next/README.md` version badge

## Acceptance

- [ ] Smoke tests pass
- [ ] No `which rp-cli` or `which codex` in any skill file
- [ ] README updated with new priority (no auto-detect)
- [ ] CHANGELOG updated
- [ ] Version bumped

## Done summary
Ran smoke tests (43/43 passed), verified detection removed from key skills. Updated README/CHANGELOG/flowctl docs to remove auto-detect from priority. Bumped version to 0.11.10. Review verdict: SHIP.
## Evidence
- Commits: 326f3e3, 09314f8
- Tests: smoke_test.sh (43/43 passed), grep verification for detection removal
- PRs: