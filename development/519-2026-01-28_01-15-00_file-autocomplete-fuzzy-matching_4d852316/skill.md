# Research – File Autocomplete Fuzzy Matching

**Date:** 2026-01-28
**Owner:** claude
**Phase:** Research

## Goal

Map the @ mention file autocomplete system and identify why fuzzy matching doesn't work (e.g., typing "pfsm" doesn't find "problem_statement_model.md").

## Findings

### Architecture Overview

```
User types "@pfsm"
    ↓
FileAutoComplete.get_search_string() → "pfsm"
    ↓
FileAutoComplete.get_candidates()
    ↓
FileFilter.complete("pfsm")  ← PROBLEM: uses prefix/substring, NOT fuzzy
    ↓
Returns [] (nothing starts with or contains "pfsm")
    ↓
AutoComplete.get_matches() → nothing to match
    ↓
Empty dropdown
```

### Relevant Files

| File | Purpose |
|------|---------|
| `src/tunacode/ui/widgets/file_autocomplete.py` | @ mention widget, extracts search string, wraps FileFilter |
| `src/tunacode/core/file_filter.py` | Core facade with default ignore patterns |
| `src/tunacode/infrastructure/file_filter.py` | Actual implementation - **uses prefix/substring, not fuzzy** |
| `.venv/.../textual_autocomplete/_autocomplete.py` | Base class, has `FuzzySearch` but only for post-filtering |
| `.venv/.../textual/fuzzy.py` | Textual's `FuzzySearch` class - available but unused |

### Root Cause

**The fuzzy matching happens in the WRONG stage of the pipeline.**

1. `FileFilter.complete(prefix)` walks the filesystem and collects files using `_matches_prefix()`
2. `_matches_prefix()` uses:
   - `startswith()` for direct children
   - `in` (substring) for nested paths
3. If nothing matches, returns `[]`
4. `AutoComplete.get_matches()` then applies fuzzy scoring to the candidates
5. But if candidates is `[]`, there's nothing to score!

**Example failure:**
```
prefix = "pfsm"
FileFilter._matches_prefix() checks:
  - "problem_statement_model.md".startswith("pfsm") → False
  - "pfsm" in "problem_statement_model.md" → False
Result: not matched, not returned
```

**What should happen:**
```
FuzzySearch.match("pfsm", "problem_statement_model.md")
  → Matches p-roblem_s-tatement_m-odel.md
  → Returns (score, [0, 8, 18, 19])
```

## Key Patterns / Solutions Found

### Option A: Use fuzzy in FileFilter (Recommended)

Replace `_matches_prefix()` logic with `textual.fuzzy.FuzzySearch`:

```python
from textual.fuzzy import FuzzySearch

class FileFilter:
    def __init__(self, ...):
        self._fuzzy = FuzzySearch(case_sensitive=False)

    def _matches_prefix(self, path: Path, name_prefix: str, search_path: Path) -> float:
        if not name_prefix:
            return 1.0  # Match all

        # For direct children, match against filename
        if path.parent == search_path:
            score, _ = self._fuzzy.match(name_prefix, path.name)
            return score

        # For nested, match against relative path
        rel_path = str(path.relative_to(search_path))
        score, _ = self._fuzzy.match(name_prefix, rel_path)
        return score
```

Then in `complete()`, collect all files with score > 0 and sort by score descending.

### Option B: Return all files, let AutoComplete filter

Change `FileFilter.complete()` to ignore `name_prefix` and return all files up to limit. Let `AutoComplete.get_matches()` do the fuzzy filtering.

**Pros:** Simpler change
**Cons:** Performance hit - walks entire tree every keystroke

### Option C: Two-pass approach

1. First pass: prefix/substring matching (fast)
2. If no results: second pass with fuzzy matching (slower but complete)

## Knowledge Gaps

- Performance implications of fuzzy matching on large codebases
- Whether `textual.fuzzy.FuzzySearch` cache is shared or per-instance
- Ideal minimum score threshold for fuzzy matches

## Implementation Notes

The fix should be in `src/tunacode/infrastructure/file_filter.py`:

1. Import `FuzzySearch` from `textual.fuzzy`
2. Create instance in `__init__`
3. Replace `_matches_prefix` to return score instead of bool
4. Update `complete()` to collect scored results and sort by score

## References

- `src/tunacode/infrastructure/file_filter.py:58-70` - `_matches_prefix()` method
- `src/tunacode/infrastructure/file_filter.py:71-127` - `complete()` method
- `.venv/.../textual/fuzzy.py:55-98` - `FuzzySearch.match()` implementation
- `.venv/.../textual_autocomplete/_autocomplete.py:511` - where fuzzy scoring happens (too late)
