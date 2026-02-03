# fn-1.3 flowctl memory commands

## Description
Add flowctl commands for manual memory management. These complement auto-capture by allowing manual additions (conventions, decisions) and retrieval.

Commands:
- `flowctl memory init` - creates `.flow/memory/` structure
- `flowctl memory add --type <type> "<content>"` - manual entry
- `flowctl memory read [--type <type>]` - dump memory
- `flowctl memory list` - show entry count per file
- `flowctl memory search "<pattern>"` - grep across memory

All commands check `memory.enabled` config first.

## Acceptance
- [ ] `flowctl memory add --type pitfall "..."` appends to pitfalls.md
- [ ] `flowctl memory add --type convention "..."` appends to conventions.md
- [ ] `flowctl memory add --type decision "..."` appends to decisions.md
- [ ] `flowctl memory read` dumps all memory files
- [ ] `flowctl memory read --type pitfalls` filters to one file
- [ ] `flowctl memory list` shows entry counts
- [ ] `flowctl memory search "pattern"` greps across files
- [ ] All commands error gracefully when memory disabled

## Done summary
- Added `flowctl memory add --type <type> "<content>"`
- Added `flowctl memory read [--type <type>]`
- Added `flowctl memory list` for entry counts
- Added `flowctl memory search "<pattern>"` for grep across files
- All commands gate on memory.enabled config

Why:
- Complements auto-capture with manual entry capability
- Enables retrieval and search of stored learnings

Verification:
- Python syntax check passed
- Manual testing of all commands passed
- flowctl validate --all passes
## Evidence
- Commits: a29de8e3c5bac283701d7bd4487d6b8b33dfc58d
- Tests: flowctl validate --all, manual testing of memory commands
- PRs: