# fn-2.2 Add flowctl codex command group

## Description

Add `flowctl codex` command group to flowctl.py, following the existing `flowctl rp` pattern.

### Commands to implement

```bash
flowctl codex check                           # Verify codex installed + version
flowctl codex impl-review <task> --base <br>  # Impl review via codex exec
flowctl codex plan-review <epic>              # Plan review via codex exec
```

### Implementation details

**`codex check`**
- Use `shutil.which("codex")` to detect availability
- Run `codex --version` to get version
- Return JSON: `{"available": true, "version": "X.Y.Z"}`

**`codex impl-review`**
- Load task info for context
- Call `gather_context_hints()` (fn-2.1) for codebase context
- Build prompt with XML structure (context_hints + review_instructions)
- **Session handling**:
  - Check for existing receipt with `session_id` → use `codex exec resume <id>`
  - No session → use `codex exec --json` and parse `thread_id` from output
- Run: `codex exec [resume <id>] --sandbox read-only --json "prompt"`
- Parse JSON output for verdict and thread_id
- Write receipt with `session_id` for re-review continuity
- Print `VERDICT=X` for eval

**`codex plan-review`**
- Read `.flow/specs/{epic_id}.md`
- Call `gather_context_hints()` for plan context
- Build prompt with XML structure (plan_spec + context_hints + review_instructions)
- **Session handling**: Same as impl-review (check receipt → resume or new)
- Run: `codex exec [resume <id>] --sandbox read-only --json "prompt"`
- Parse JSON, extract verdict, write receipt with session_id

**Session continuity design**:
- First review: `codex exec --json` → parse `thread_id` from `{"type":"thread.started","thread_id":"..."}`
- Store in receipt: `{"session_id": "019baa19-...", ...}`
- Re-review: Read receipt → `codex exec resume <session_id> "prompt"`
- Fallback: If resume fails (session expired), start new session
- **Never use `--last`** - conflicts with parallel usage

### Key patterns to follow

- `/Users/gordon/work/gmickel-claude-marketplace/plugins/flow-next/scripts/flowctl.py:157-173` - `require_rp_cli()` pattern
- `/Users/gordon/work/gmickel-claude-marketplace/plugins/flow-next/scripts/flowctl.py:2660-2726` - `cmd_rp_setup_review` atomic pattern
- `/Users/gordon/work/gmickel-claude-marketplace/plugins/flow-next/scripts/flowctl.py:3109-3185` - argparse setup

### Files to modify

- `plugins/flow-next/scripts/flowctl.py` - add codex command group
## Acceptance
- [ ] `flowctl codex check --json` returns `{"available": true/false, "version": "X.Y.Z"}`
- [ ] `flowctl codex impl-review fn-1.1 --base main --json` executes codex review
- [ ] `flowctl codex plan-review fn-1 --json` executes codex plan review
- [ ] Verdict extracted from output: SHIP or NEEDS_WORK
- [ ] Receipt includes `"mode": "codex"` and `"session_id": "..."`
- [ ] Re-review uses `codex exec resume <session_id>` when receipt exists
- [ ] Fallback to new session if resume fails
- [ ] Graceful error when codex not installed
- [ ] Existing `flowctl rp` commands still work (no regression)
## Done summary
- Added `flowctl codex` command group with check, impl-review, plan-review
- Session continuity via thread_id storage/resume
- Helper functions: require_codex, get_codex_version, run_codex_exec, parse_codex_thread_id, parse_codex_verdict, build_review_prompt

Why:
- Enables cross-platform code reviews using Codex CLI
- Alternative to RepoPrompt for Windows/Linux users

Verification:
- smoke_test.sh passes
- `flowctl codex check --json` returns {"available": true, "version": "0.80.0"}
- Help output correct for all subcommands
- Helper functions tested: verdict/thread_id parsing, prompt building
## Evidence
- Commits: a06f103cab5cbe2950b950cea71a7a3e5302ed14
- Tests: smoke_test.sh, codex check --json
- PRs: