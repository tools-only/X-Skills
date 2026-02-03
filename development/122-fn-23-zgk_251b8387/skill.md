# Readable Epic IDs (Slugified Titles)

## Overview

Change epic ID format from random suffix (`fn-5-x7k`) to slugified title (`fn-5-typed-validators`). The numeric prefix stays for ordering/uniqueness; only the suffix changes from random to descriptive.

**Current**: `fn-23-zgk` (random 3-char suffix)
**New**: `fn-23-readable-epic-ids` (slugified title)

## Scope

**In scope:**
- Add `slugify()` function to flowctl.py (stdlib only, Django pattern)
- Update `generate_epic_suffix()` to use slugified title
- Update regex in `parse_id()` to allow variable-length slugs
- Update 30+ error messages with new format examples
- Update CLI help text
- Update all skill/doc files with new ID examples

**Out of scope:**
- Migration of existing IDs (they remain valid via backwards compat)
- Renaming existing epics
- Collision resolution (number prefix ensures uniqueness)

## Approach

### Slugify implementation (Django pattern, no deps)

```python
import re
import unicodedata

def slugify(text: str, max_length: int = 40) -> str:
    """Convert text to URL-safe slug. Max 40 chars for ID suffix."""
    text = str(text)
    text = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode("ascii")
    text = re.sub(r"[^\w\s-]", "", text.lower())
    text = re.sub(r"[-\s]+", "-", text).strip("-_")
    if max_length and len(text) > max_length:
        text = text[:max_length].rsplit("-", 1)[0]  # truncate at word boundary
    return text or None  # return None if empty for fallback
```

### ID generation change

```python
# Current (flowctl.py:2664-2665)
suffix = generate_epic_suffix()  # random 3-char
epic_id = f"fn-{epic_num}-{suffix}"

# New
slug = slugify(args.title, max_length=40)
suffix = slug if slug else generate_epic_suffix()  # fallback to random if empty
epic_id = f"fn-{epic_num}-{suffix}"
```

### Regex update (flowctl.py:590)

```python
# Current
r"^fn-(\d+)(?:-[a-z0-9]{3})?(?:\.(\d+))?$"

# New (supports: fn-N, fn-N-xxx, fn-N-longer-slug)
r"^fn-(\d+)(?:-[a-z0-9][a-z0-9-]*[a-z0-9]|-[a-z0-9]{1,3})?(?:\.(\d+))?$"
```

## Quick commands

```bash
# Test slugify
python3 -c "
import re, unicodedata
def slugify(t, m=40):
    t = unicodedata.normalize('NFKD', str(t)).encode('ascii', 'ignore').decode('ascii')
    t = re.sub(r'[^\w\s-]', '', t.lower())
    t = re.sub(r'[-\s]+', '-', t).strip('-_')
    return t[:m].rsplit('-', 1)[0] if m and len(t) > m else t
print(slugify('Add OAuth Login'))  # add-oauth-login
print(slugify('Fix bug #123'))     # fix-bug-123
print(slugify(''))                 # (empty)
"

# Test regex
python3 -c "
import re
p = r'^fn-(\d+)(?:-[a-z0-9][a-z0-9-]*[a-z0-9]|-[a-z0-9]{1,3})?(?:\.(\d+))?\$'
for t in ['fn-1', 'fn-1-abc', 'fn-1-add-oauth', 'fn-1.2', 'fn-1-abc.2', 'fn-1-add-oauth.2']:
    print(f'{t}: {bool(re.match(p, t))}')"

# Run smoke tests
cd /tmp && /path/to/plugins/flow-next/scripts/smoke_test.sh
```

## Acceptance

- [ ] `flowctl epic create --title "Add OAuth"` produces `fn-N-add-oauth`
- [ ] `flowctl epic create --title ""` falls back to random suffix
- [ ] `flowctl epic create --title "!@#$"` falls back to random suffix
- [ ] Existing `fn-N` and `fn-N-xxx` IDs still parse correctly
- [ ] `fn-N-longer-slug` format parses correctly
- [ ] All 45 smoke tests pass
- [ ] Error messages show new format examples
- [ ] CLI help text shows new format
- [ ] All skill docs updated with new ID examples

## References

- Django slugify: `django/utils/text.py` (stdlib only)
- Current regex: `flowctl.py:590`
- Suffix generation: `flowctl.py:579-582`
- Epic creation: `flowctl.py:2664-2665`
- Error messages: 30+ locations (grep "Expected format")
