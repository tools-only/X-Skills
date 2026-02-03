# fn-13-pxj.6 Update documentation

## Description
Update documentation to cover new async control commands.

**NOTE**: Only update existing files. Do NOT create new documentation files.

### Files to Update

**1. `plugins/flow-next/README.md`**

Add "Controlling Ralph" section:
- `flowctl status` usage and output examples
- `flowctl ralph pause/resume/stop/status` commands
- Sentinel file mechanism (for manual control)
- `flowctl task reset` for rollback/retry workflows
- Example: pause → fix requirement → resume workflow

Also update **Troubleshooting** section:
- Replace outdated `flowctl task set --status pending` reference with `flowctl task reset`
- The current docs mention `task set` which doesn't exist in flowctl.py

**2. Root `CLAUDE.md`**

No changes needed - it references README which will be updated.

### Documentation Style
- Keep examples concise
- Show both CLI and manual (file-based) approaches where applicable
- Include JSON output examples for automation use cases
- Ensure all documented commands actually exist in flowctl.py

### Key Files
- `plugins/flow-next/README.md`
## Acceptance
- [ ] README.md has "Controlling Ralph" section
- [ ] Examples show `flowctl status`, `flowctl ralph`, `flowctl task reset`
- [ ] Troubleshooting section updated (replace `task set` with `task reset`)
- [ ] All documented commands exist in flowctl.py
- [ ] Documentation is concise and matches actual command behavior
- [ ] No new documentation files created
## Done summary
- Added "Controlling Ralph" subsection with CLI + sentinel file docs
- Fixed Troubleshooting: replaced nonexistent `task set` with `task reset`
- Removed jq dependency from example code

Why: Document new async control features, fix broken documentation

Verification: CI 40/40, README renders correctly
## Evidence
- Commits: 7b9969e9557088074a5c108c0814d98c182ca67a
- Tests: plugins/flow-next/scripts/ci_test.sh
- PRs: