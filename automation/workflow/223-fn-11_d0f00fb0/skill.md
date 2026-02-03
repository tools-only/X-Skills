# fn-1.1 Config + memory directory structure

## Description
Add `.flow/config.json` support to flowctl for project settings. Create config commands and memory init that creates `.flow/memory/` structure when enabled.

Key implementation:
- `load_flow_config()` helper that handles missing file gracefully
- `get_config(key, default)` for nested keys like `memory.enabled`
- Both `flowctl init` and `/flow-next:setup` create config if missing
- Memory commands check `memory.enabled` before running

## Acceptance
- [ ] `flowctl config get <key>` returns nested config values
- [ ] `flowctl config set <key> <value>` updates config.json
- [ ] Default config has `memory.enabled: false`
- [ ] `flowctl init` creates config.json with defaults
- [ ] `/flow-next:setup` creates config.json with defaults (idempotent)
- [ ] Setup Step 7 mentions: "Memory: disabled. Enable with flowctl config set memory.enabled true"
- [ ] `flowctl memory init` creates `.flow/memory/` (only if enabled)
- [ ] Memory templates created: pitfalls.md, conventions.md, decisions.md
- [ ] Memory commands error gracefully when disabled

## Done summary
- Added `load_flow_config()` and `get_config(key, default)` helpers
- Implemented `flowctl config get/set` commands for nested keys
- Implemented `flowctl memory init` with gating on memory.enabled
- Updated `flowctl init` to create config.json with defaults
- Updated setup workflow Step 7 to mention memory system

Why:
- Provides foundation for memory system opt-in behavior
- All memory operations gated behind single config flag

Verification:
- flowctl validate --all passes
- Manual testing of config get/set commands
- Manual testing of memory init (enabled/disabled states)
## Evidence
- Commits: 4398a30be610c53384bf07d3e646d82f771f70ec
- Tests: flowctl validate --all
- PRs: