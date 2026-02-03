# fn-13-pxj.8 Add flowctl epic add-dep/rm-dep commands

## Description
Add `flowctl epic add-dep` and `flowctl epic rm-dep` commands to manage epic-level dependencies via CLI.

Currently `depends_on_epics` can only be set by editing JSON directly. The existing `dep add` command only works for tasks.

### Files to Change

**1. `plugins/flow-next/scripts/flowctl.py`**

Add subcommands under `epic`:
```python
# Around L1550-1600 where epic subparsers are defined
p_epic_add_dep = epic_sub.add_parser("add-dep", help="Add epic dependency")
p_epic_add_dep.add_argument("epic", help="Epic ID")
p_epic_add_dep.add_argument("depends_on", help="Epic ID to depend on")
p_epic_add_dep.add_argument("--json", action="store_true")
p_epic_add_dep.set_defaults(func=cmd_epic_add_dep)

p_epic_rm_dep = epic_sub.add_parser("rm-dep", help="Remove epic dependency")
p_epic_rm_dep.add_argument("epic", help="Epic ID")
p_epic_rm_dep.add_argument("depends_on", help="Epic ID to remove from deps")
p_epic_rm_dep.add_argument("--json", action="store_true")
p_epic_rm_dep.set_defaults(func=cmd_epic_rm_dep)
```

Add handlers:
```python
def cmd_epic_add_dep(args):
    epic_data = load_epic(args.epic)
    dep_data = load_epic(args.depends_on)  # Validate exists
    
    if args.depends_on == args.epic:
        error_exit("Epic cannot depend on itself", args.json)
    
    deps = epic_data.get("depends_on_epics", [])
    if args.depends_on in deps:
        # Already exists, no-op success
        json_output({"success": True, ...}) if args.json else print(...)
        return
    
    deps.append(args.depends_on)
    epic_data["depends_on_epics"] = deps
    epic_data["updated_at"] = now_iso()
    save_epic(epic_data)
    # Output

def cmd_epic_rm_dep(args):
    epic_data = load_epic(args.epic)
    deps = epic_data.get("depends_on_epics", [])
    
    if args.depends_on not in deps:
        # Not in deps, no-op success
        json_output({"success": True, ...}) if args.json else print(...)
        return
    
    deps.remove(args.depends_on)
    epic_data["depends_on_epics"] = deps
    epic_data["updated_at"] = now_iso()
    save_epic(epic_data)
    # Output
```

**2. `plugins/flow-next/README.md`**

Add to CLI reference section (~L535):
```markdown
| `flowctl epic add-dep <epic> <dep>` | Add epic-level dependency |
| `flowctl epic rm-dep <epic> <dep>` | Remove epic-level dependency |
```

**3. `plugins/flow-next/docs/flowctl.md`**

Add command reference (~L50):
```markdown
## flowctl epic add-dep
Add an epic-level dependency.

flowctl epic add-dep <epic-id> <depends-on-epic-id> [--json]

Example:
flowctl epic add-dep fn-2 fn-1 --json
{"success": true, "id": "fn-2", "depends_on_epics": ["fn-1"]}

## flowctl epic rm-dep
Remove an epic-level dependency.

flowctl epic rm-dep <epic-id> <depends-on-epic-id> [--json]
```

**4. `plugins/flow-next/scripts/smoke_test.sh`**

Update `depends_on_epics gate` test (~L584-611) to use CLI:
```bash
# Replace Python JSON edit with:
scripts/flowctl epic add-dep "$DEP_CHILD_ID" "$DEP_BASE_ID" --json
```

**5. `plugins/flow-next/scripts/ci_test.sh`**

Add test:
```bash
test_epic_add_dep() {
    epic1=$($FLOWCTL epic create --title "Base epic" --json | python3 -c '...')
    epic2=$($FLOWCTL epic create --title "Child epic" --json | python3 -c '...')
    
    $FLOWCTL epic add-dep "$epic2" "$epic1" --json
    deps=$($FLOWCTL show "$epic2" --json | python3 -c 'import json,sys; print(",".join(json.load(sys.stdin)["depends_on_epics"]))')
    [[ "$deps" == "$epic1" ]] && pass "epic add-dep" || fail "epic add-dep"
    
    $FLOWCTL epic rm-dep "$epic2" "$epic1" --json
    deps=$($FLOWCTL show "$epic2" --json | python3 -c 'import json,sys; print(",".join(json.load(sys.stdin)["depends_on_epics"]))')
    [[ -z "$deps" ]] && pass "epic rm-dep" || fail "epic rm-dep"
}
```

### Key Files
- `plugins/flow-next/scripts/flowctl.py` - add commands (~L1550, ~L1700)
- `plugins/flow-next/README.md` - CLI reference table
- `plugins/flow-next/docs/flowctl.md` - command docs
- `plugins/flow-next/scripts/smoke_test.sh` - update test
- `plugins/flow-next/scripts/ci_test.sh` - add test
## Acceptance
- [ ] `flowctl epic add-dep <epic> <dep>` adds dependency
- [ ] `flowctl epic rm-dep <epic> <dep>` removes dependency
- [ ] Self-dependency rejected with error
- [ ] Non-existent epic rejected with error
- [ ] Duplicate add is no-op (success)
- [ ] Remove non-existent dep is no-op (success)
- [ ] `--json` flag works for both commands
- [ ] README CLI reference updated
- [ ] docs/flowctl.md command reference added
- [ ] smoke_test.sh updated to use CLI
- [ ] ci_test.sh test added and passes
## Done summary
- Added `cmd_epic_add_dep()` to add epic-level dependencies
- Added `cmd_epic_rm_dep()` to remove epic-level dependencies
- Validation: self-dep error, non-existent epic error
- Idempotent: duplicate add/remove non-existent is no-op success

Why: Enable CLI management of epic deps (currently requires JSON editing)

Verification: Manual tests for add/rm/errors; CI 31/31
Follow-up: CI test + doc updates in tasks .5/.6
## Evidence
- Commits: 47f74c62247d40acc2bfe23688f81dfbda25150b
- Tests: plugins/flow-next/scripts/ci_test.sh
- PRs: