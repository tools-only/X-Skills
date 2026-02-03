# fn-17-wti.1 Add review.backend to flowctl default config

## Description
Add `review.backend` to the default config structure in `flowctl.py`.

## File to modify

`plugins/flow-next/scripts/flowctl.py`

## Changes

1. Find `get_default_config()` function (around line 76-78)
2. Add `review` key with `backend` subkey:

```python
def get_default_config() -> dict:
    """Return default config structure."""
    return {
        "memory": {"enabled": False},
        "planSync": {"enabled": False},
        "review": {"backend": None}
    }
```

3. The `None` value means "not configured yet" - skills will error if this is None and no flag/env provided.

## Why None instead of "auto"

- Explicit is better than implicit
- Forces user to make a choice via setup or flag
- Prevents LLM deviation from "auto-detect" logic

## Test

```bash
# After change, init a fresh .flow/ and check config
rm -rf /tmp/test-config && mkdir /tmp/test-config && cd /tmp/test-config
python3 /path/to/flowctl.py init --json
cat .flow/config.json | jq '.review'
# Should output: {"backend": null}
```
## Acceptance
- [ ] `get_default_config()` returns dict with `review.backend` key
- [ ] Value is `None` (Python None, becomes JSON null)
- [ ] Existing configs upgraded via `flowctl init` get the new key merged in
- [ ] Smoke tests pass
## Done summary
Added review.backend key with None value to get_default_config() in flowctl.py. Existing configs get the new key via deep_merge on upgrade.
## Evidence
- Commits: 3a3248c906b97a0ef66dd9ed3b3a41b487f73925
- Tests: flowctl init --json (fresh dir), config get review.backend (upgrade test), smoke_test.sh (43 tests pass)
- PRs: