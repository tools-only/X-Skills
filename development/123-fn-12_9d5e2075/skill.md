# fn-1.2 Capture hook

## Description
Extend ralph-guard.py PostToolUse handler to capture learnings from NEEDS_WORK reviews. When chat-send returns NEEDS_WORK/MAJOR_RETHINK, parse feedback and append to pitfalls.md.

Key implementation:
- Check `memory.enabled` config before processing
- `extract_feedback(response)` parses review into structured format
- `is_learnable(feedback)` filters to actionable patterns only
- Append to `.flow/memory/pitfalls.md` with date, task, issue, fix, category

Filter criteria (`is_learnable`):
- Has specific actionable fix (not vague)
- References code pattern, API, or convention
- Not a one-off typo or obvious bug

## Acceptance
- [ ] PostToolUse handler detects chat-send with NEEDS_WORK
- [ ] Checks `memory.enabled` before processing
- [ ] Parses feedback into structured format (issue, fix, category)
- [ ] `is_learnable()` filters out non-actionable items
- [ ] Appends to pitfalls.md with correct format
- [ ] Handles missing memory dir gracefully

## Done summary
- Added memory capture in ralph-guard.py PostToolUse handler
- extract_feedback() parses NEEDS_WORK reviews into structured issues
- is_learnable() filters to actionable patterns (framework, API, convention)
- classify_issue() categorizes entries for memory format
- append_to_pitfalls() writes entries to .flow/memory/pitfalls.md

Why:
- Captures high-signal learnings from reviewer feedback
- Only activates when memory.enabled is true

Verification:
- Python syntax check passed
- Unit test of extract_feedback and is_learnable passed
- flowctl validate --all passes
## Evidence
- Commits: 1471f8ca7bb4373bd576e057a8aefc1283e22ff3
- Tests: python3 -m py_compile, unit test of memory functions
- PRs: