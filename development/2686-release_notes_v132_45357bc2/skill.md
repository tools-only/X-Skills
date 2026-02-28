# Lár v1.3.2 Release Notes

**"The Developer Experience Update"**

This release addresses the documentation gap from v1.3.0 by adding comprehensive usage examples and documentation for `AuditLogger` and `TokenTracker`.

## What's New

### Documentation Improvements

#### 1. New Example: Custom Logger/Tracker
**File:** `examples/patterns/16_custom_logger_tracker.py`

Demonstrates:
- Default automatic creation (recommended for most use cases)
- Custom dependency injection for advanced workflows
- Shared tracker for cost aggregation across multiple executors
- Direct access to logger history and token summaries

**Usage:**
```bash
python examples/patterns/16_custom_logger_tracker.py
```

#### 2. Enhanced v1.3.0 Release Notes
Added complete "How to Use" section with:
- Option 1: Automatic (default behavior)
- Option 2: Custom injection (advanced)
- Code examples for both patterns
- Use cases for custom injection

#### 3. Updated README.md
- Added observability features section to `GraphExecutor` documentation
- Included code examples for default and advanced usage
- Added reference to new example in metacognition features list

#### 4. Updated Documentation
**File:** `docs/core-concepts/3-audit-log.md`

New section: "Modular Observability (v1.3.0)"
- Detailed `AuditLogger` usage patterns
- Detailed `TokenTracker` usage patterns
- Real-world examples for cost aggregation
- GxP-compliant audit trail examples

## Impact

**Before v1.3.2:**
- Users had to read source code to understand `AuditLogger` and `TokenTracker`
- No documented examples of custom injection
- Confusion about when to use custom instances vs defaults

**After v1.3.2:**
- ✅ Clear usage examples in release notes, README, and docs
- ✅ Dedicated example file with 4 usage patterns
- ✅ Guidance on when to use each approach
- ✅ Significantly improved developer experience

## No Breaking Changes

This is purely a documentation release. No API changes, no code changes to the core framework.

All existing code continues to work exactly as before.

## Upgrade

```bash
pip install --upgrade lar-engine
```

Or with Poetry:
```bash
poetry add lar-engine@^1.3.2
```

## Documentation Coverage

| Location | Status |
|----------|--------|
| Release Notes (v1.3.0) | ✅ Updated with usage examples |
| Examples | ✅ New 16_custom_logger_tracker.py |
| README.md | ✅ GraphExecutor section updated |
| Docs | ✅ core-concepts/3-audit-log.md updated |

---

**Developer Experience:** Significantly improved for v1.3.0 observability features.
