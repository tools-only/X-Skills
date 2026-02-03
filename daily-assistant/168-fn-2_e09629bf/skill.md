# Codex Review Backend

## Overview

Add Codex CLI as alternative review backend for cross-model code reviews. Enables Windows/Linux users (who can't use RepoPrompt) to get multi-model review quality gates.

## Scope

- New `flowctl codex` command group (check, impl-review, plan-review)
- Context hints via inline implementation (repo-scout pattern, not subagent)
- Update skills to branch on backend (impl-review, plan-review)
- Update Ralph templates (config.env, prompt_plan.md, prompt_work.md)
- Update ralph-guard.py hooks for codex patterns
- Add smoke tests and e2e test variant

## Approach

### Key Design Decisions

1. **Use `codex exec` for all reviews** - `codex review --base` doesn't accept custom prompts (Issue #7825)
2. **Context hints inline in flowctl** - can't spawn subagents from codex, use grep-based approach
3. **Session continuity via `codex exec resume`** - store session_id in receipt, resume on re-review
4. **Config priority**: env var > .flow/config.json > config.env default

### Session Continuity (Re-review Loop)

RP maintains chat context automatically (same window). Codex needs explicit session tracking:

1. **First review**: `codex exec --json ...` → parse `thread_id` from output
2. **Store in receipt**: `{"session_id": "019baa19-...", ...}`
3. **Re-review**: Read receipt → `codex exec resume <session_id> "prompt"`

This gives RP-like context - reviewer sees previous feedback + changes.

**Never use `--last`** - conflicts with parallel codex usage or multiple projects.

### Command Structure

```bash
flowctl codex check                           # Verify codex installed + auth
flowctl codex impl-review <task> --base <br>  # Impl review via codex exec
flowctl codex plan-review <epic>              # Plan review via codex exec
```

### Prompt Structure

```
<context_hints>
Consider these related files:
- path/to/file.ts:42 - reason
</context_hints>

<review_instructions>
[Carmack prompt]
Output <verdict>SHIP</verdict> or <verdict>NEEDS_WORK</verdict>
</review_instructions>
```

## Quick Commands

```bash
# Run smoke tests
cd /tmp && /Users/gordon/work/gmickel-claude-marketplace/plugins/flow-next/scripts/smoke_test.sh
cd /tmp && /Users/gordon/work/gmickel-claude-marketplace/plugins/flow-next/scripts/ralph_smoke_test.sh

# Test codex check (if codex installed)
FLOWCTL="plugins/flow-next/scripts/flowctl"
$FLOWCTL codex check --json
```

## Acceptance

- [ ] `flowctl codex check` returns availability + version
- [ ] `flowctl codex impl-review` sends review to codex, extracts verdict, writes receipt
- [ ] `flowctl codex plan-review` sends plan review, extracts verdict, writes receipt
- [ ] Session continuity works (re-review uses `codex exec resume`)
- [ ] Skills branch correctly based on backend config
- [ ] Ralph works with `PLAN_REVIEW=codex WORK_REVIEW=codex`
- [ ] ralph-guard.py validates codex calls (blocks direct codex without wrapper, blocks --last)
- [ ] Smoke tests pass with codex backend
- [ ] Graceful fallback when codex unavailable
- [ ] Docs updated: READMEs, ralph-e2e-notes.md, ralph-getting-started.md
- [ ] RepoPrompt still recommended as primary, Codex as cross-platform alternative

## Risks

1. **Codex auth**: Non-interactive auth check limited to `codex --version`
2. **Different model behavior**: GPT may interpret Carmack prompt differently than Claude
3. **Session expiry**: Codex sessions may expire - need fallback to new session with context summary

## References

- Spec: plans/codex-review-backend.md
- Existing RP implementation: flowctl.py:2474-2726
- repo-scout pattern: agents/repo-scout.md
- Ralph hooks: scripts/hooks/ralph-guard.py
