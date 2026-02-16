# COV-010: Test worker state edge cases

## Overview
Create tests for get_worker_state with invalid IDs, set_worker_state with None in tests/unit/test_state_workers.py.

## Files to Create
- `tests/unit/test_state_workers.py`

## Files to Modify
None

## Files to Read
- `zerg/state.py`

## Requirements

Implement the task as described. Follow TDD approach:
1. Write failing tests first
2. Implement code to make tests pass
3. Refactor if needed

## Verification

```bash
pytest tests/unit/test_state_workers.py -v --tb=short
```

## Acceptance Criteria
- [ ] Test get_worker_state with non-existent worker ID
- [ ] Test set_worker_state with None fields
- [ ] Test worker state update preserves existing fields
- [ ] Test get_all_workers returns empty when no workers
- [ ] Test worker state serialization roundtrip
