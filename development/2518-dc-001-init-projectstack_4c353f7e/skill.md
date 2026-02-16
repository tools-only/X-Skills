# DC-001: Update init.py to use ProjectStack

**Level**: 1 | **Critical Path**: No | **Estimate**: 20 min

## Objective

Replace the single-type `detect_project_type()` function with `ProjectStack` from `security_rules.py` to support multi-language project detection.

## Files Owned

- `zerg/commands/init.py` (modify)

## Files to Read

- `zerg/security_rules.py` (reference ProjectStack, detect_project_stack)

## Implementation Steps

1. Import `detect_project_stack` and `ProjectStack` from `zerg.security_rules`
2. Replace `detect_project_type()` call with `detect_project_stack(Path("."))`
3. Update `create_config()` to accept `ProjectStack` instead of `str | None`
4. Update console output to show all detected languages, frameworks, databases
5. Update `show_summary()` to display multi-language info
6. Maintain backwards compatibility with config structure

## Code Changes

```python
# Old
project_type = detect_project_type()
if project_type:
    console.print(f"Detected project type: [cyan]{project_type}[/cyan]")

# New
from zerg.security_rules import detect_project_stack, ProjectStack

stack = detect_project_stack(Path("."))
if stack.languages:
    langs = ", ".join(sorted(stack.languages))
    console.print(f"Detected languages: [cyan]{langs}[/cyan]")
if stack.frameworks:
    frameworks = ", ".join(sorted(stack.frameworks))
    console.print(f"Detected frameworks: [cyan]{frameworks}[/cyan]")
```

## Verification

```bash
cd /tmp && rm -rf test-dc && mkdir test-dc && cd test-dc
touch requirements.txt package.json go.mod
zerg init --no-security-rules 2>&1 | grep -E 'python|javascript|go'
# Should show all three languages
```

## Acceptance Criteria

- [ ] `detect_project_stack()` replaces `detect_project_type()`
- [ ] Init output shows all detected languages
- [ ] Init output shows detected frameworks if any
- [ ] Config still saves project_type for backwards compat (use primary language)
- [ ] No ruff errors: `ruff check zerg/commands/init.py`
