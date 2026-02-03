# fn-12 Collision-resistant epic IDs (fn-N-xxx format)

## Overview
Change epic ID format from `fn-N` to `fn-N-xxx` where `xxx` is a 3-char alphanumeric suffix. Prevents merge conflicts when multiple team members create epics simultaneously (both getting `fn-5` → collision).

## Problem
Sequential numeric IDs cause:
- Merge conflicts when team members create epics in parallel
- ID collisions requiring manual renumbering
- File conflicts in `.flow/epics/`, `.flow/tasks/`, `.flow/specs/`

## Solution
Format: `fn-N-xxx` where:
- `N` = next available number (scan-based, same as today)
- `xxx` = 3-char random alphanumeric (a-z0-9), e.g., `fn-5-x7k`

Task IDs follow: `fn-N-xxx.M` (e.g., `fn-5-x7k.1`)

## Approach
**Dual-pattern support** - no migration needed:
- Accept both `fn-N` (legacy) and `fn-N-xxx` (new) formats
- New epics always get suffix
- Existing epics continue working unchanged

## Files to modify

### Core (flowctl.py)
- `parse_id()` line ~336: Add suffix group to regex
- `scan_max_epic_id()` line ~933: Extract number ignoring suffix
- `epic create` line ~1494: Generate suffix on creation
- `EpicId`/`TaskId` types: Update validation

### Ralph hooks (ralph-guard.py)
- Receipt parsing ~line 50: Update `plan-` and `impl-` patterns
- ID extraction: Handle suffix in receipt filenames

### TUI (flow-next-tui)
- `runs.ts`: Update `TASK_ID_PATTERN` and epic extraction regex
- Display: Truncate suffix in tight spaces if needed

### Tests
- `test_flowctl.py`: Update hardcoded IDs and expected patterns
- `test_runs.ts`: Update TUI test fixtures

### Documentation
- `plugins/flow-next/README.md`: Document new format
- `CHANGELOG.md`: Note backwards-compatible change
- `plans/` runbooks: Update examples

## Quick commands
- `bun test` (TUI tests)
- `plugins/flow-next/scripts/smoke_test.sh`
- `plugins/flow-next/scripts/ralph_smoke_test.sh`

## Acceptance
- [ ] `flowctl epic create "test"` produces `fn-N-xxx` format
- [ ] `flowctl show fn-5` still works (legacy)
- [ ] `flowctl show fn-5-x7k` works (new)
- [ ] Ralph receipts work with new format
- [ ] TUI displays and filters new format correctly
- [ ] All smoke tests pass
- [ ] No migration required for existing `.flow/` dirs

## Risk assessment

**High confidence (95%+):**
- Dual-pattern regex is straightforward: `fn-(\d+)(?:-[a-z0-9]{3})?`
- Existing `fn-N` IDs will match (optional suffix group)
- No migration = no data risk
- 28 smoke + 15 ralph tests catch regressions

**Medium confidence areas:**
- `scan_max_epic_id()` must extract number ignoring suffix - easy but needs care
- TUI regex updates - isolated, has tests

**Low risk but needs attention:**
- Receipt filename parsing in ralph-guard.py
- Ensuring glob patterns (`fn-*.json`) still work (they will)

**What could go wrong:**
1. Typo in regex breaks existing ID parsing → caught by smoke tests
2. scan_max returns wrong number → caught by epic create test
3. TUI can't parse new format → caught by bun test

**Verdict:** High confidence. Dual-pattern approach is inherently safe. Every change has test coverage.

## References
- Analysis: Deep search of all ID patterns in codebase
- Patterns found: 10+ regex across flowctl.py, ralph-guard.py, TUI
- Backwards compat: Dual-pattern matching, no breaking changes
