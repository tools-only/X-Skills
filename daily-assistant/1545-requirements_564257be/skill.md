# Requirements: Security Hardening

**Status**: APPROVED
**Created**: 2026-02-04
**Issues**: #131, #135, #142

## Problem Statement

ZERG has three security vulnerabilities related to command injection and path traversal that need remediation:

1. **Command Injection via shell=True** (Issue #131, Severity: HIGH, CWE-78)
   - Three modules bypass `CommandExecutor` and use `subprocess.run(..., shell=True)` directly
   - Affected: `HypothesisTestRunner`, `StepExecutor`, `RecoveryPlanner`

2. **Weak Command Prefix Validation** (Issue #135, Severity: MEDIUM, CWE-183)
   - `CommandExecutor.validate_command()` uses simple prefix matching
   - `python -c` and `npx` in allowlist enable arbitrary code execution
   - `trust_commands=True` flag conditionally enables shell=True

3. **Path Traversal in Security Scan** (Issue #142, Severity: MEDIUM, CWE-22)
   - `run_security_scan()` follows symlinks via `os.walk()` without validation
   - Can scan files outside intended directory boundary

## Functional Requirements

### FR-1: Migrate shell=True Callers to CommandExecutor

**FR-1.1**: Refactor `HypothesisTestRunner.test()` (hypothesis_engine.py:194-200)
- Replace `subprocess.run(shell=True)` with `CommandExecutor.execute()`
- Remove SAFE_PREFIXES list (CommandExecutor handles validation)
- Preserve timeout and output capture behavior

**FR-1.2**: Refactor `StepExecutor._execute_step()` (step_executor.py:272-278)
- Use `CommandExecutor` for all command execution
- Initialize executor with `trust_commands=True` for task-graph.json verification commands
- Add validation before execution

**FR-1.3**: Refactor `RecoveryPlanner.execute_step()` (recovery.py:367-374)
- Replace `subprocess.run(shell=True)` with `CommandExecutor.execute()`
- Apply `shlex.quote()` to template variables before interpolation
- Preserve timeout (SUBPROCESS_TIMEOUT=10) behavior

### FR-2: Harden CommandExecutor Allowlist

**FR-2.1**: Remove dangerous prefixes from default allowlist
- Remove: `python -c`, `python3 -c`, `npx`
- These enable arbitrary code execution and bypass the safety intent

**FR-2.2**: Add pattern-based validation for specific use cases
- Create `SAFE_PYTHON_PATTERNS` for permitted python -c patterns
- Allow only specific, validated patterns (e.g., JSON operations, pytest)

**FR-2.3**: Deprecate `trust_commands` flag
- Log warning when `trust_commands=True` is used
- Document migration path for existing callers
- Add TODO for eventual removal

### FR-3: Prevent Path Traversal in Security Scan

**FR-3.1**: Disable symlink following in `os.walk()`
- Change to `os.walk(scan_path, followlinks=False)`
- This is explicit defense even though False is the default

**FR-3.2**: Validate resolved paths stay within scan boundary
- After resolving each file path, check `resolved.is_relative_to(scan_path)`
- Skip and log warning for symlinks pointing outside boundary

**FR-3.3**: Handle symlink resolution errors gracefully
- Catch OSError, ValueError during path resolution
- Skip problematic paths with logging

## Non-Functional Requirements

### NFR-1: Backward Compatibility
- All existing tests must pass
- Task-graph.json verification commands must continue to work
- Recovery plan commands must execute successfully

### NFR-2: Security Audit Trail
- All command executions logged via CommandExecutor audit_log
- Symlink violations logged with warning level

### NFR-3: Performance
- No measurable performance regression
- Path validation adds O(1) per file overhead

### NFR-4: Test Coverage
- Add unit tests for each fixed vulnerability
- Add integration tests for CommandExecutor usage in new callers
- Add test for symlink boundary enforcement

## Out of Scope

- Redesigning CommandExecutor from scratch
- Adding sandboxing or containerization for command execution
- Implementing full AppArmor/seccomp profiles
- Changes to Docker security configuration

## Acceptance Criteria

1. **AC-1**: `grep -r "shell=True" zerg/` returns only `command_executor.py` (conditional path)
2. **AC-2**: `CommandExecutor` is imported in hypothesis_engine.py, step_executor.py, recovery.py
3. **AC-3**: `python -c` and `npx` removed from `ALLOWED_COMMAND_PREFIXES`
4. **AC-4**: `os.walk(..., followlinks=False)` used in security.py
5. **AC-5**: All existing tests pass: `pytest tests/`
6. **AC-6**: New tests added for each vulnerability fix
7. **AC-7**: bandit scan passes: `bandit -r zerg/ -ll`

## Verification Commands

```bash
# Check shell=True usage is consolidated
grep -rn "shell=True" zerg/ | grep -v "command_executor.py"
# Expected: empty output (only command_executor.py should have conditional shell)

# Check dangerous prefixes removed
grep -E "(python -c|python3 -c|npx)" zerg/command_executor.py
# Expected: empty output

# Check followlinks=False
grep -n "os.walk" zerg/security.py
# Expected: contains followlinks=False

# Run tests
pytest tests/ -v

# Run security linter
bandit -r zerg/ -ll
```

## Technical Notes

### Template Variable Safety (FR-1.3)
Recovery templates use `{feature}` and `{worker_id}` placeholders:
```python
# Before (vulnerable to injection via feature name)
cmd = f"zerg logs --worker {worker_id}"

# After (safe)
import shlex
cmd = f"zerg logs --worker {shlex.quote(worker_id)}"
```

### CommandExecutor Integration Pattern
```python
from zerg.command_executor import CommandExecutor, get_executor

# Option 1: Get default executor
executor = get_executor()
result = executor.execute(command)

# Option 2: Create with specific config
executor = CommandExecutor(
    working_dir=working_dir,
    timeout=30,
    trust_commands=False  # Default, safest
)
```

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Regression in task verification | Medium | High | Extensive test coverage |
| Breaking recovery commands | Low | Medium | Test recovery plan execution |
| Performance degradation | Low | Low | Benchmark before/after |

## Dependencies

- No external dependencies
- Uses existing `zerg.command_executor` module
- Uses standard library `shlex`, `pathlib`

## Open Questions

None - requirements are derived from detailed GitHub issue specifications.
