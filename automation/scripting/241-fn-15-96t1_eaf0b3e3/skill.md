# fn-15-96t.1 Add planSync.enabled config to flowctl.py

## Description
TBD

## Acceptance
- [ ] TBD

## Done summary
Added planSync.enabled config key to flowctl.py get_default_config(). Defaults to false (opt-in). Verified set/get works.
## Evidence
- Commits: 42dc44f533454e7277e7e62d19e466c2dab284c2
- Tests: flowctl config set planSync.enabled true --json, flowctl config get planSync.enabled --json
- PRs: