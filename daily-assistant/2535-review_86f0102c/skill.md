# ZERG Review

Three-stage code review workflow (Spec, Quality, Security).

## Usage

```bash
/zerg:review [--mode prepare|self|receive|full] [--no-security]
```

## Modes

### prepare
Prepare PR for review:
- Generate change summary
- Check spec compliance
- Create review checklist

### self
Self-review checklist:
- Code compiles without errors
- All tests pass locally
- No hardcoded values or secrets
- Error handling is appropriate
- Edge cases are handled
- Code is readable and well-named

### receive
Process review feedback:
- Parse review comments
- Track addressed items
- Generate response

### full (default)
Complete three-stage review:
1. Spec compliance check
2. Code quality review
3. Security scan (via consolidated `run_security_scan()` engine)

## Stage 3: Security

Runs the consolidated security engine (`run_security_scan()`) against modified/reviewed files. Checks 15 capability areas including secrets detection, injection patterns, crypto misuse, CVE dependency scanning, authentication flaws, access control issues, SSRF patterns, deserialization risks, path traversal, XSS, error handling, logging gaps, hardcoded credentials, insecure configuration, and sensitive data exposure.

### --no-security flag

Skip the security scan stage. Use with caution — prints a WARNING when invoked. Useful for quick spec/quality-only reviews where security has already been verified separately via `/zerg:security`.

## Examples

```bash
# Full two-stage review
/zerg:review

# Prepare for PR
/zerg:review --mode prepare

# Self-review checklist
/zerg:review --mode self
```

## Output

```
Code Review Results
========================================
Status: PASSED
Files Reviewed: 5

Stage 1 (Spec):     ✓
Stage 2 (Quality):  ✓
Stage 3 (Security): ✓  (15 capabilities, 0 findings)
```

## Task Tracking

On invocation, create a Claude Code Task to track this command:

Call TaskCreate:
  - subject: "[Review] Code review {mode}"
  - description: "Running code review in {mode} mode for current feature."
  - activeForm: "Running code review"

Immediately call TaskUpdate:
  - taskId: (the Claude Task ID)
  - status: "in_progress"

On completion, call TaskUpdate:
  - taskId: (the Claude Task ID)
  - status: "completed"

## Help

When `--help` is passed in `$ARGUMENTS`, display usage and exit:

```
/zerg:review — Three-stage code review workflow (Spec → Quality → Security).

Flags:
  --mode MODE           Review mode: prepare|self|receive|full (default: full)
  --no-security         Skip Stage 3 security scan (prints WARNING)
  --help                Show this help message
```

## Exit Codes

- 0: Review passed
- 1: Issues found
- 2: Review error
