# Issue Standards

This document outlines the preferred conventions for creating and maintaining GitHub issues in this project. These are personal preferences of the maintainer, not strict requirements. Following these guidelines helps keep issues focused and discussions productive.

## Title Format

Issues should use one of these prefixes:

| Type | Prefix | Example |
|------|--------|---------|
| Bug | `[Bug]:` | `[Bug]: Update command hangs on uv-only installations` |
| Feature | `[Feature]:` | `[Feature]: Model picker lacks recent models history` |
| Question | `[Question]:` | `[Question]: Pydantic removal contract discussion` |

### Title Guidelines

- Describe the **problem**, not the solution
- Use present tense ("lacks", "hangs", "fails")
- Be specific but concise

**Good:** `[Bug]: Tool call JSON parsing silently fails without model feedback`
**Bad:** `[Bug]: Fix JSON parsing to add error feedback`

## Issue Body Structure

The issue body should describe the **problem only**. Solutions belong in follow-up comments.

### Recommended Sections

```markdown
## Summary
One paragraph describing what is broken or missing.

## Context
Why this matters. What triggers this issue. When users encounter it.

## Root Cause (for bugs)
Technical explanation of why this happens.
```

### Optional Sections

```markdown
## Impact
What users experience. What breaks.

## Related Issues
Links to related issues if applicable.
```

### Better as Comments (Not in Body)

These are better placed in follow-up comments rather than the issue body:

- Implementation code
- "Proposed Solution" sections
- "Action Items" or task checklists
- "Files to Modify" lists
- Acceptance criteria

## Solutions Go in Comments

After creating the issue, add a **comment** with:

- Proposed solution
- Implementation approach
- Files to modify
- Acceptance criteria
- Code examples

This keeps the issue focused on the problem and allows multiple solution proposals.

### Example Flow

1. Create issue with problem description
2. Add comment: "Proposed approach: ..."
3. Discussion happens in comments
4. Implementation PR references the issue

## Labels

Apply appropriate labels after creating the issue:

| Label | When to Use |
|-------|-------------|
| `bug` | Something is broken |
| `enhancement` | New feature or improvement |
| `documentation` | Docs only changes |
| `good first issue` | Simple, well-scoped tasks |
| `help wanted` | Community contributions welcome |
| `hard` | Complex or risky changes |

## Examples

### Good Issue

```markdown
Title: [Bug]: Update command hangs on uv-only installations

## Summary
The `/update` flow appears to hang because the update check shells out to
`pip index` with no timeout and no uv-aware fallback.

## Context
Triggered when users run `/update` or `/update check`. Also affects
background update check at UI shutdown.

## Root Cause
`check_for_updates()` uses `subprocess.run(["pip", "index", ...])` without
a timeout. On uv-only installs, pip may not exist or behave differently.

## Impact
Users see the UI stall. uv-only installs may misreport "Already on latest".
```

Then in a **comment**:
```markdown
## Proposed Fix
1. Add 5s timeout to subprocess call
2. Add uv fallback: `uv pip show tunacode`
3. Use semver comparison instead of string comparison

## Files
- `src/tunacode/utils/system/paths.py:161`
- `src/tunacode/ui/commands/__init__.py:425`
```

### Bad Issue

```markdown
Title: [Feature]: Add timeout to update command

## Summary
We should add a timeout and uv support.

## Implementation
1. Change line 161 to add timeout=5
2. Add try/except for uv
3. Install packaging library for semver

## Files to Modify
- paths.py
- commands/__init__.py
```

This could be improved:
- Title describes solution, not problem
- Body jumps straight to implementation
- Missing context about why this matters
