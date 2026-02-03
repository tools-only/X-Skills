# fn-17-wti.4 Remove runtime detection from flow-next-plan-review

## Description
Remove runtime `which rp-cli` / `which codex` detection from plan-review skill files.

## Files to modify

1. `plugins/flow-next/skills/flow-next-plan-review/SKILL.md`
2. `plugins/flow-next/skills/flow-next-plan-review/workflow.md`

## Changes to SKILL.md

### Remove detection block (lines ~42-51)

Delete:
```bash
# Check available backends
HAVE_RP=$(which rp-cli >/dev/null 2>&1 && echo 1 || echo 0)
HAVE_CODEX=$(which codex >/dev/null 2>&1 && echo 1 || echo 0)

# Get configured backend
BACKEND="${FLOW_REVIEW_BACKEND:-}"
if [[ -z "$BACKEND" ]]; then
  BACKEND="$($FLOWCTL config get review.backend 2>/dev/null | jq -r '.value // empty')"
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

### Remove "If no backend configured and both available" section (lines ~53-68)

Delete the interactive question that asks user to choose backend.
Replace with nothing (handled by error below).

### Remove fallback detection (lines ~72-79)

Delete:
```bash
# Fallback to available
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

### Update Backend Selection docs (lines ~20-27)

Update priority list from:
```
1. `--review=rp|codex|export|none` argument
2. `FLOW_REVIEW_BACKEND` env var (`rp`, `codex`, `none`)
3. `.flow/config.json` → `review.backend`
4. Interactive prompt if both rp-cli and codex available (and not in Ralph mode)
5. Default: whichever is available (rp preferred)
```

To:
```
1. `--review=rp|codex|export|none` argument
2. `FLOW_REVIEW_BACKEND` env var (`rp`, `codex`, `none`)
3. `.flow/config.json` → `review.backend`
4. **ERROR** - "No review backend configured"
```

## Changes to workflow.md

### Remove Phase 0 detection (lines ~20-38)

Delete:
```bash
# Check available backends
HAVE_RP=$(which rp-cli >/dev/null 2>&1 && echo 1 || echo 0)
HAVE_CODEX=$(which codex >/dev/null 2>&1 && echo 1 || echo 0)

# Get configured backend (priority: env > config)
BACKEND="${FLOW_REVIEW_BACKEND:-}"
if [[ -z "$BACKEND" ]]; then
  BACKEND="$($FLOWCTL config get review.backend 2>/dev/null | jq -r '.value // empty' 2>/dev/null || echo "")"
fi

# Fallback to available (rp preferred)
if [[ -z "$BACKEND" ]]; then
  if [[ "$HAVE_RP" == "1" ]]; then BACKEND="rp"
  elif [[ "$HAVE_CODEX" == "1" ]]; then BACKEND="codex"
  else BACKEND="none"; fi
fi
```

Replace with:
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

## Done summary
Removed runtime `which rp-cli` / `which codex` detection from plan-review skill files (SKILL.md and workflow.md). Now uses config-only priority: --review flag > env > config > error. Error message guides users to /flow-next:setup or --review flag.
## Evidence
- Commits: a3f73d520510e68b61a720422c66e327d3b4b991
- Tests: grep -r 'which rp-cli' plugins/flow-next/skills/flow-next-plan-review/ (no matches)
- PRs: