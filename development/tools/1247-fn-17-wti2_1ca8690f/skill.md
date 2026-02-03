# fn-17-wti.2 Add review backend question to /flow-next:setup

## Description
Add a review backend question to the `/flow-next:setup` workflow.

## File to modify

`plugins/flow-next/skills/flow-next-setup/workflow.md`

## Changes

### Step 6: Add detection before questions

Before the `AskUserQuestion` block, add tool detection:

```bash
# Detect available review backends
HAVE_RP=$(which rp-cli >/dev/null 2>&1 && echo 1 || echo 0)
HAVE_CODEX=$(which codex >/dev/null 2>&1 && echo 1 || echo 0)

# Read current config value if exists
CURRENT_BACKEND=$("${PLUGIN_ROOT}/scripts/flowctl" config get review.backend --json 2>/dev/null | jq -r '.value // empty')
```

### Step 6: Add question to AskUserQuestion block

Add a new question to the existing `questions` array (after Plan-Sync, before Docs):

```json
{
  "header": "Review",
  "question": "Which review backend for Carmack-level reviews?",
  "options": [
    {"label": "Codex CLI (Recommended)", "description": "Cross-platform, uses GPT 5.2 High. Detected: <yes/no>"},
    {"label": "RepoPrompt", "description": "macOS only, visual context builder. Detected: <yes/no>"},
    {"label": "None", "description": "Skip reviews, can configure later with --review flag"}
  ],
  "multiSelect": false
}
```

**Dynamic options based on detection:**
- If `HAVE_CODEX=0`: Change Codex description to include "(not detected)"
- If `HAVE_RP=0`: Change RepoPrompt description to include "(not detected)"
- If neither detected: Add note "Neither rp-cli nor codex detected. Install one for review support."

**Pre-selection for re-runs:**
- If `CURRENT_BACKEND` is set, note it in the question: "Current: <backend>"

### Step 7: Process answer

Add handling for Review answer:

```bash
# Review backend
case "$review_answer" in
  "Codex"*) REVIEW_BACKEND="codex" ;;
  "RepoPrompt"*) REVIEW_BACKEND="rp" ;;
  *) REVIEW_BACKEND="none" ;;
esac

"${PLUGIN_ROOT}/scripts/flowctl" config set review.backend "$REVIEW_BACKEND" --json
```

### Step 8: Update summary

Add to the Configuration section:

```
Configuration:
- Memory: <enabled|disabled>
- Plan-Sync: <enabled|disabled>
- Review backend: <codex|rp|none>
```
## Acceptance
- [ ] Setup workflow detects rp-cli and codex availability
- [ ] Review backend question added to AskUserQuestion block
- [ ] Question shows detection status (detected/not detected)
- [ ] Re-run setup shows current config value
- [ ] Answer is written to config via `flowctl config set review.backend`
- [ ] Summary shows configured review backend
- [ ] Works when neither tool detected (defaults to "none")
## Done summary
Added review backend question to /flow-next:setup workflow with detection for rp-cli and codex, answer mapping to config, and summary output. Includes jq dependency note.
## Evidence
- Commits: 97f49822c8b4b4f8dc7d3f0baec4b9e6cd58eb0e, f2cb6f59e0c9cd1b1c5cf7fe152946837321d8f2
- Tests: Manual inspection of workflow.md
- PRs: