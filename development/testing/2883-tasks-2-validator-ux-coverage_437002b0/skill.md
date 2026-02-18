---
feature: plugin-validator UX and coverage gaps
status: complete
created: 2026-02-14
source: .claude/BACKLOG.md (P1)
target-file: ./plugins/plugin-creator/scripts/plugin_validator.py
test-dir: ./plugins/plugin-creator/tests/
pre-existing-failures: 18
---

# Tasks: plugin-validator UX and coverage gaps

## Feature Summary

Fix UX bugs and expand validation coverage in `plugin_validator.py`.
Four independent sub-issues from groomed backlog item.

## Pre-existing Test State

**18 failed, 260 passed, 1 skipped** as of 2026-02-14.
Do NOT regress any passing tests. Do NOT fix pre-existing failures unless they are directly in the scope of a task below.

## Context Manifest

### Primary Files

- `./plugins/plugin-creator/scripts/plugin_validator.py` (3045 lines)
- `./plugins/plugin-creator/tests/conftest.py` (341 lines)

### Architecture

- **Validator Protocol** (line 249): `validate(path) -> ValidationResult`, `can_fix() -> bool`, `fix(path) -> list[str]`
- **FileType enum** (line 138): `SKILL`, `AGENT`, `COMMAND`, `PLUGIN`, `UNKNOWN`
- **Error codes** (lines 66-112): FM001-FM010, SK001-SK007, LK001-LK002, PD001-PD003, PL001-PL005, NR001-NR002
- **Validator selection** (lines 2885-2925): Skills/Agents/Commands get FrontmatterValidator, NameFormatValidator, DescriptionValidator, NamespaceReferenceValidator. Skills additionally get ComplexityValidator, InternalLinkValidator, ProgressiveDisclosureValidator. Plugins get PluginStructureValidator.
- **Result collection** (lines 2928-2931): Iterates validators (not files), creating duplicate path entries
- **Report display** (lines 2620-2710): Prints per-result, counts results as files

### Test Conventions

- Module loaded via `importlib.util.spec_from_file_location("plugin_validator", ...)` (conftest.py)
- Fixtures: `cli_runner`, `sample_skill_dir`, `sample_agent_dir`, `sample_plugin_dir`, `mock_frontmatter_file`
- Pattern: create tmp files with frontmatter, instantiate validator, call `.validate()`, assert on `ValidationResult`
- Class-based test organization, `pytest.mark.parametrize`

### Linting

- Ruff py311 strict mode
- Per-file ignores for `plugins/plugin-creator/scripts/*.py`: B008, TC003
- Per-file ignores for `**/scripts/**`: ANN401, DOC, PLC0415, PLR0911, S, T201
- Per-file ignores for `**/tests/**`: ANN, D, DOC, E501, EXE, N, PLC, PLR, S, SLF, T

---

## Task 1: Fix report counting — unique files not validator invocations

**Status**: ✅ COMPLETE
**Started**: 2026-02-14T16:21:00Z
**Completed**: 2026-02-14T16:30:00Z
**Priority**: P0 (Must Have)
**Estimated scope**: ~50 lines changed

### Problem

Lines 2928-2931 iterate over `validators` (not files), appending `(path, result)` for each validator. When 1 file has 4 validators, `results` list has 4 entries pointing to the same path. The report at line 2630 prints "PASSED" 4 times for 1 file. The summary at line 2702 shows "Total files: 4" for 1 file.

### Changes Required

1. Restructure result collection (lines 2928-2931) to group results by file path
2. Update report display (lines 2620-2710) to:
   - Show each file once with all validator results beneath it
   - Label each validator result with the validator class name (e.g., "FrontmatterValidator: PASSED")
   - A file passes overall only if ALL its validators pass
3. Update summary (line 2702) to count unique file paths

### Acceptance Criteria

- Running validator on 1 file shows "Total files: 1"
- Each validator result is labeled with the validator name
- Summary counts unique files, not validator invocations
- All currently-passing tests remain passing

### Test Requirements

- Add test: single file with multiple validators shows "Total files: 1"
- Add test: validator names appear in per-file output

---

## Task 2: File-type-aware DescriptionValidator — SK004/SK005 scoping

**Status**: ✅ COMPLETE
**Started**: 2026-02-14T16:31:00Z
**Completed**: 2026-02-14T16:40:00Z
**Priority**: P0 (Must Have)
**Estimated scope**: ~80 lines changed

### Problem

`DescriptionValidator.validate()` (line 1790) applies SK004 ("description too short") and SK005 ("missing trigger phrases") to all file types. Commands have a different frontmatter schema and don't need trigger phrases.

### Changes Required

1. Add `file_type: FileType` parameter to `DescriptionValidator.__init__()` (or pass FileType to `validate()`)
2. SK005 ("missing trigger phrases"): Only fire on SKILL files
3. SK004 ("description too short"): Fire on SKILL and AGENT files, not COMMAND
4. Update validator selection (lines 2885-2925) to pass FileType when constructing DescriptionValidator
5. Define new error code series CM001+ for command-specific description checks (future — stub only)

### Acceptance Criteria

- SK005 only fires on SKILL files (not COMMAND or AGENT)
- SK004 fires on SKILL and AGENT files but not COMMAND
- Validator selection passes FileType context to DescriptionValidator
- All currently-passing tests remain passing

### Test Requirements

- Add test: command file does NOT receive SK005 warning
- Add test: command file does NOT receive SK004 warning
- Add test: agent file does NOT receive SK005 warning
- Add test: agent file still receives SK004 warning
- Add test: skill file still receives both SK004 and SK005

---

## Task 3: Hook validation — FileType.HOOK and HookValidator

**Status**: ✅ COMPLETE
**Started**: 2026-02-14T16:41:00Z
**Completed**: 2026-02-14T16:55:00Z
**Priority**: P1 (Should Have)
**Estimated scope**: ~200 lines new code

### Problem

`FileType` enum (line 138) has no HOOK variant. `detect_file_type()` returns UNKNOWN for `.js` files in `hooks/` directories and `hooks.json` files.

### Changes Required

1. Add `HOOK_SCRIPT` and `HOOK_CONFIG` to `FileType` enum (line 138)
2. Update `detect_file_type()` to recognize:
   - `.js` files in `hooks/` directories -> HOOK_SCRIPT
   - `hooks.json` files -> HOOK_CONFIG
3. Define new error code constants: HK001 (invalid hooks.json), HK002 (invalid event type), HK003 (invalid hook structure)
4. Create `HookValidator` class implementing Validator protocol:
   - For HOOK_CONFIG: validate JSON structure, valid event types, valid matcher patterns
   - For HOOK_SCRIPT: validate file is valid JavaScript (basic syntax check), has proper shebang
5. Update validator selection to assign HookValidator for HOOK_SCRIPT and HOOK_CONFIG files
6. Reference: `/plugin-creator:claude-hooks-reference-2026` skill for hooks.json schema

### Acceptance Criteria

- `detect_file_type()` recognizes `.js` files in `hooks/` directories
- `detect_file_type()` recognizes `hooks.json` files
- `HookValidator` validates hooks.json structure
- New HK001-HK003 error codes defined and used
- All currently-passing tests remain passing

### Test Requirements

- New test file: `test_hook_validator.py`
- Test: hooks.json with valid structure passes
- Test: hooks.json with invalid JSON fails (HK001)
- Test: hooks.json with invalid event type fails (HK002)
- Test: .js hook file detected as HOOK_SCRIPT
- Test: hooks.json detected as HOOK_CONFIG
- Test: detect_file_type for hook files returns correct FileType

---

## Task 4: Remove dead code — nested skill error message references

**Status**: ✅ COMPLETE
**Started**: 2026-02-14T16:15:00Z
**Completed**: 2026-02-14T16:20:00Z
**Priority**: P2 (Could Have)
**Estimated scope**: ~10 lines changed

### Problem

**CORRECTION from backlog**: Lines 904-911 are NOT dead code. The nested pattern `plugins/{plugin}/skills/*/{name}/SKILL.md` exists in the repo (e.g., `python3-development/skills/testing/*/SKILL.md`, `python3-development/skills/development/*/SKILL.md`).

However, the error message strings at lines 758-760 and 770-773 reference the nested pattern with a glob wildcard (`skills/*/`) which is confusing in user-facing error messages. These should be rephrased for clarity.

### Changes Required

1. Update error message at lines 758-760 to clarify the nested pattern (e.g., "plugins/{plugin}/skills/{name}/SKILL.md or plugins/{plugin}/skills/{category}/{name}/SKILL.md")
2. Update error message at lines 770-773 similarly
3. Verify no test depends on exact error message text (grep tests for the old strings)

### Acceptance Criteria

- Error messages clearly describe both flat and nested skill paths
- No wildcard globs in error messages shown to users
- All currently-passing tests remain passing

### Test Requirements

- Verify existing NR tests still pass after message change
- No new tests needed (cosmetic change)

---

## Dependencies

Tasks 1-4 are independent and can be implemented in parallel.
Task 3 is the largest (new validator class, new error codes, new test file).

## Execution Order (suggested)

1. Task 4 (smallest, cosmetic) — builds confidence, fast win
2. Task 1 (UX fix) — most visible improvement
3. Task 2 (file-type awareness) — false positive elimination
4. Task 3 (hook validation) — largest, new feature
