# fn-13-pxj.7 Bump version to 0.8.0 and update CHANGELOG

## Description
Bump version to 0.8.0 and update CHANGELOG.

### Files to Update

**1. Version files (manual, no bump script for flow-next):**
- `plugins/flow-next/.claude-plugin/plugin.json` → version: "0.8.0"
- `.claude-plugin/marketplace.json` → flow-next version in plugins array
- `plugins/flow-next/README.md` → Version badge
- `README.md` → Flow-next badge

**2. CHANGELOG.md:**
Add `[flow-next 0.8.0]` entry with:
- `flowctl status` command for monitoring
- `flowctl ralph pause/resume/stop/status` for run control
- `flowctl task reset [--cascade]` for retry workflows
- Sentinel file mechanism in ralph.sh
- CI tests for all new commands

### Validation
```bash
jq . .claude-plugin/marketplace.json
jq . plugins/flow-next/.claude-plugin/plugin.json
```

### Key Files
- `plugins/flow-next/.claude-plugin/plugin.json`
- `.claude-plugin/marketplace.json`
- `plugins/flow-next/README.md`
- `README.md`
- `CHANGELOG.md`
## Acceptance
- [ ] plugin.json version is "0.8.0"
- [ ] marketplace.json flow-next version is "0.8.0"
- [ ] README badges updated to 0.8.0
- [ ] CHANGELOG.md has [flow-next 0.8.0] entry
- [ ] Entry lists all new commands
- [ ] JSON files validate with jq
## Done summary
- Bumped version to 0.8.0 in plugin.json, marketplace.json, README badges
- Added CHANGELOG entry with all new features

Why: Release new version with async control features

Verification: JSON validates with jq, CI 40/40
## Evidence
- Commits: 466b0c49cf569a23c90828d7d41be793199cec51
- Tests: plugins/flow-next/scripts/ci_test.sh
- PRs: