# fn-17-wti.5 Remove runtime detection from flow-next-work

## Description
Remove runtime `which rp-cli` / `which codex` detection from flow-next-work skill.

## File to modify

`plugins/flow-next/skills/flow-next-work/SKILL.md`

## Changes

### Remove detection block (lines ~58-68)

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
# Check configured backend (priority: env > config)
CONFIGURED_BACKEND="${FLOW_REVIEW_BACKEND:-}"
if [[ -z "$CONFIGURED_BACKEND" ]]; then
  CONFIGURED_BACKEND="$($FLOWCTL config get review.backend 2>/dev/null | jq -r '.value // empty')"
fi
```

### Update question logic (lines ~96-167)

**Key change**: Questions no longer depend on `HAVE_RP` or `HAVE_CODEX`. Instead:

1. **If CONFIGURED_BACKEND is set** (from env or config): Only ask branch question
2. **If CONFIGURED_BACKEND is empty**: Ask both branch AND review questions, but the review question offers all options without detection status

Replace the complex "If both rp-cli AND codex available" / "If only rp-cli available" / "If only codex available" / "If neither available" logic with:

### If backend already configured (env or config)

Skip review question, only ask branch (but show override hint):
```
Quick setup: Where to work?
a) Current branch  b) New branch  c) Isolated worktree

(Reply: "a", "current", or just tell me)
(Tip: --review=rp|codex|none overrides configured backend)
```

### If no backend configured

Ask both questions:
```
Quick setup before starting:

1. **Branch** — Where to work?
   a) Current branch
   b) New branch
   c) Isolated worktree

2. **Review** — Run Carmack-level review after?
   a) Codex CLI
   b) RepoPrompt
   c) Export for external LLM
   d) None (configure later with --review flag)

(Reply: "1a 2a", "current branch, codex", or just tell me naturally)
```

**If user picks review option but tool not installed**: Error at review time (not here).

### Update defaults (lines ~160-167)

Change from:
```
- Review = configured backend if set, else `codex` if available, else `rp` if available, else `none`
```

To:
```
- Review = configured backend if set, else error at review time
```

## Acceptance

- [ ] No `which rp-cli` or `which codex` in SKILL.md
- [ ] No `HAVE_RP` or `HAVE_CODEX` variables
- [ ] Review question appears when no backend configured
- [ ] Review question skipped when backend configured (env or config)
- [ ] Default behavior: if no review choice and no config, skip review (user can pass --review flag)
- [ ] Override hint shown when backend configured: "(Tip: --review=rp|codex|none overrides configured backend)"

## Done summary
Removed runtime `which rp-cli`/`which codex` detection from flow-next-work skill. Simplified question logic: ask review question only when no backend configured (env or config), show override hint otherwise.
## Evidence
- Commits: cbc75d45a1af57e72614b24a2197f5ce32b26e16
- Tests: grep verification: no which rp-cli/codex or HAVE_* in SKILL.md
- PRs: