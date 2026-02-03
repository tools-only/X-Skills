# fn-1.6 Documentation

## Description
Document memory system in relevant places. Key point: memory is a flow-next feature (not Ralph-specific), works in manual + Ralph modes.

Clarify:
- Config in `.flow/config.json` (not Ralph's `env.config`)
- Enable with `flowctl config set memory.enabled true`
- Auto-capture only happens via Ralph hook
- Manual additions via `flowctl memory add`

## Acceptance
- [ ] Memory section added to flow-next README
- [ ] Memory mentioned in ralph.md (note: separate config)
- [ ] CLAUDE.md updated with memory quick-reference
- [ ] Clear distinction: flow-next feature vs Ralph-only
- [ ] flowctl memory commands documented with examples

## Done summary
- Added Memory System section to flow-next README
- Updated .flow/ directory tree to show config.json and memory files
- Added Memory Capture section to ralph.md
- Added memory quick-reference to CLAUDE.md

Key clarifications:
- Memory is flow-next feature, not Ralph-specific
- Config in .flow/config.json (separate from Ralph's config.env)
- Works in manual + Ralph modes

Verification:
- flowctl validate --all passes
## Evidence
- Commits: aea13e84163a72c8fba7380df747720d335046a0
- Tests: flowctl validate --all
- PRs: