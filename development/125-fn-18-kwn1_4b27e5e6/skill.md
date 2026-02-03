# fn-18-kwn.1 Update flowctl rp to auto-open RP window via rp-cli open

## Description

`rp-cli open <path>` can now launch RepoPrompt windows programmatically. Update `flowctl rp pick-window` to use this when no existing window is found.

## Implementation

1. **Update `pick_window()` in flowctl.py**
   ```python
   def pick_window(repo_root: str) -> int:
       # Try to find existing window
       result = subprocess.run(['rp-cli', 'list-windows', '--json'], capture_output=True)
       if result.returncode == 0:
           windows = json.loads(result.stdout)
           for w in windows:
               if w.get('path') == repo_root:
                   return w['id']

       # No window found - open one
       result = subprocess.run(['rp-cli', 'open', repo_root, '--json'], capture_output=True)
       if result.returncode == 0:
           data = json.loads(result.stdout)
           return data['window_id']

       raise FlowctlError(f"Failed to open RepoPrompt: {result.stderr}")
   ```

2. **Update error messaging**
   - Remove "please open RepoPrompt first" errors
   - Add helpful error if rp-cli not installed

3. **Test with Ralph**
   - Verify RP review flow works without pre-opened window

## Key Files

- `plugins/flow-next/scripts/flowctl.py` - rp subcommand implementation

## Acceptance

- [ ] `flowctl rp pick-window --repo-root .` returns window ID (opens if needed)
- [ ] Reuses existing window if repo already open
- [ ] Clear error if rp-cli not installed
- [ ] Ralph RP reviews work without manual window setup

## Done summary
TBD

## Evidence
- Commits:
- Tests:
- PRs:
