# fn-17-wti.3 Remove runtime detection from flow-next-impl-review

## Description
Remove runtime `which rp-cli` / `which codex` detection from impl-review skill files.

## Files to modify

1. `plugins/flow-next/skills/flow-next-impl-review/SKILL.md`
2. `plugins/flow-next/skills/flow-next-impl-review/workflow.md`

## Changes to SKILL.md

### Remove detection block (lines ~43-52)

Delete:
```bash
HAVE_RP=$(which rp-cli >/dev/null 2>&1 && echo 1 || echo 0)
HAVE_CODEX=$(which codex >/dev/null 2>&1 && echo 1 || echo 0)

# Check configured backend (priority: env > config)
CONFIGURED_BACKEND="${FLOW_REVIEW_BACKEND:-}"
if [[ -z "$CONFIGURED_BACKEND" ]]; then
  CONFIGURED_BACKEND="$($FLOWCTL config get review.backend 2>/dev/null | jq -r '.value // empty')"
fi
```

Replace with simpler config-only read:
```bash
# Priority: --review flag > env > config
BACKEND="${FLOW_REVIEW_BACKEND:-}"
if [[ -z "$BACKEND" ]]; then
  BACKEND="$($FLOWCTL config get review.backend 2>/dev/null | jq -r '.value // empty')"
fi
```

### Remove fallback detection (lines ~75-78)

Delete:
```bash
if [[ -z "$BACKEND" ]]; then
  if [[ "$HAVE_RP" == "1" ]]; then BACKEND="rp"
  elif [[ "$HAVE_CODEX" == "1" ]]; then BACKEND="codex"
  else BACKEND="none"; fi
fi
```

Replace with error:
```bash
if [[ -z "$BACKEND" || "$BACKEND" == "null" ]]; then
  echo "Error: No review backend configured."
  echo "Run /flow-next:setup to configure, or pass --review=rp|codex|none"
  exit 1
fi
```

### Update "If no backend configured" section (lines ~53-68)

Remove the interactive question that asks user to choose backend.
Replace with: error message telling user to run setup or pass flag.

## Changes to workflow.md

### Remove Phase 0 detection (lines ~21-35)

Delete detection and fallback, replace with:
```bash
# Priority: --review flag > env > config (flag parsed in SKILL.md)
BACKEND="${FLOW_REVIEW_BACKEND:-}"
if [[ -z "$BACKEND" ]]; then
  BACKEND="$($FLOWCTL config get review.backend 2>/dev/null | jq -r '.value // empty' 2>/dev/null || echo "")"
fi

if [[ -z "$BACKEND" || "$BACKEND" == "null" ]]; then
  echo "Error: No review backend configured."
  echo "Run /flow-next:setup to configure, or pass --review=rp|codex|none"
  exit 1
fi

echo "Review backend: $BACKEND (override: --review=rp|codex|none)"
```
## Acceptance
- [ ] No `which rp-cli` or `which codex` in SKILL.md
- [ ] No `which rp-cli` or `which codex` in workflow.md
- [ ] No `HAVE_RP` or `HAVE_CODEX` variables
- [ ] Priority order: --review flag > FLOW_REVIEW_BACKEND env > config
- [ ] Clear error message when no backend configured
- [ ] Error message mentions both /flow-next:setup and --review flag
- [ ] Success output shows override hint: "Review backend: X (override: --review=rp|codex|none)"
- [ ] Smoke tests pass
## Done summary
Removed runtime `which rp-cli` and `which codex` detection from impl-review skill. Backend selection now uses flag > env > config priority with clear error when unconfigured.
## Evidence
- Commits: 4fde7a6f6d822f5c4528cd0f0fe0853b496a9c6b
- Tests: plugins/flow-next/scripts/smoke_test.sh
- PRs: