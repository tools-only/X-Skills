# fn-17-wti.6 Remove runtime detection from flow-next-plan

## Description
Remove runtime `which rp-cli` / `which codex` detection from flow-next-plan skill.

## File to modify

`plugins/flow-next/skills/flow-next-plan/SKILL.md`

## Changes

### Remove detection block (lines ~43-53)

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

### Update question logic (lines ~69-141)

**Key change**: Questions no longer depend on `HAVE_RP` or `HAVE_CODEX`. Instead:

1. **If CONFIGURED_BACKEND is set** (from env or config): Only ask research question (if relevant)
2. **If CONFIGURED_BACKEND is empty**: Ask both research AND review questions

Replace the complex "If both rp-cli AND codex available" / "If only rp-cli available" / "If only codex available" logic with:

### If backend already configured (env or config)

Skip review question, only ask research (but show override hint):
```
Quick setup: Use RepoPrompt for deeper context?
a) Yes, context-scout (slower, thorough)
b) No, repo-scout (faster)

(Reply: "a", "b", or just tell me)
(Tip: --review=rp|codex|none overrides configured backend)
```

### If no backend configured

Ask both questions:
```
Quick setup before planning:

1. **Research approach** — Use RepoPrompt for deeper context?
   a) Yes, context-scout (slower, thorough)
   b) No, repo-scout (faster)

2. **Review** — Run Carmack-level review after?
   a) Codex CLI
   b) RepoPrompt
   c) Export for external LLM
   d) None (configure later with --review flag)

(Reply: "1a 2a", "1b 2d", or just tell me naturally)
```

**Special case**: Research question only makes sense if rp-cli might be used. Since we no longer detect tools, always show both options for research. Context-scout will error at runtime if rp-cli not available.

### Update defaults (lines ~132-141)

Change from:
```
- Review = configured backend if set, else `codex` if available, else `rp` if available, else `none`
```

To:
```
- Review = configured backend if set, else `none` (user should run setup or pass --review flag)
```

## Acceptance

- [ ] No `which rp-cli` or `which codex` in SKILL.md
- [ ] No `HAVE_RP` or `HAVE_CODEX` variables
- [ ] Review question appears when no backend configured
- [ ] Review question skipped when backend configured (env or config)
- [ ] Research question always shows both options (context-scout errors at runtime if rp-cli missing)
- [ ] Override hint shown when backend configured: "(Tip: --review=rp|codex|none overrides configured backend)"

## Done summary
Removed runtime `which rp-cli` / `which codex` detection from flow-next-plan skill. Simplified question logic to use configured backend from env/config, with override hint when backend is pre-configured.
## Evidence
- Commits: d960c3ff2e27a77fbc6f5ce5e3f54bf1f27f5c07, 38fb2022def6bf548c72436ac4a327b56990231a
- Tests: grep detection check
- PRs: