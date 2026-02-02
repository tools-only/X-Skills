# Conversion Examples

Before/after examples for each resource type.

## Example 1: Code Reviewer Subagent

### Before (`.claude/agents/go-reviewer.md`)

```yaml
---
name: go-reviewer
description: Reviews Go code for correctness, patterns, and security
model: claude-opus-4-20250514
allowed-tools:
  - Read
  - Grep
  - Glob
readonly: true
---

# Go Code Reviewer Agent

You are an expert Go code reviewer. Analyze code for:
1. Correctness and potential bugs
2. Go idioms and best practices
3. Security vulnerabilities
4. Performance issues

Provide specific, actionable feedback with code examples.
```

### After (`.claude/skills/go-reviewer/SKILL.md`)

```yaml
---
name: go-reviewer
description: |
  Reviews Go code for correctness, patterns, and security.
  Use after writing or modifying Go files, or when reviewing PRs.
  Identifies bugs, security issues, and Go idiom violations.
context: fork
model: opus
allowed-tools: Read Grep Glob
---

# Go Code Reviewer

You are an expert Go code reviewer. Analyze code for:
1. Correctness and potential bugs
2. Go idioms and best practices
3. Security vulnerabilities
4. Performance issues

Provide specific, actionable feedback with code examples.

## Review Process

1. Read the file(s) to review
2. Check for common Go issues
3. Verify error handling patterns
4. Assess security implications
5. Suggest improvements with examples
```

**Key changes:**
- Added `context: fork` for subagent isolation
- Simplified `model: claude-opus-4-20250514` → `opus`
- Enhanced description with trigger keywords
- Converted `readonly: true` to explicit `allowed-tools: Read Grep Glob`

---

## Example 2: Test Runner Command

### Before (`.claude/commands/run-tests.md`)

```yaml
---
description: Run project tests with coverage
argument-hint: <test-pattern>
---

Run tests matching the pattern: $ARGUMENTS

1. Execute: `go test -v -cover ./... -run "$1"`
2. Report coverage summary
3. List any failing tests with details
```

### After (`.claude/skills/run-tests/SKILL.md`)

```yaml
---
name: run-tests
description: Run project tests with coverage. Invoke with /run-tests <pattern>.
disable-model-invocation: true
argument-hint: <test-pattern>
---

# Run Tests

Execute tests matching the provided pattern with coverage reporting.

## Usage

Invoke with: `/run-tests <pattern>`

## Process

1. Execute: `go test -v -cover ./... -run "$1"`
2. Report coverage summary
3. List any failing tests with details

## Arguments

- `$ARGUMENTS` - Full argument string
- `$1` - Test pattern to match
```

**Key changes:**
- Added `name` field (required)
- Added `disable-model-invocation: true` (user-invoked only)
- Documented argument variables

---

## Example 3: Always-Apply Style Rule

### Before (`.cursor/rules/code-style.mdc`)

```yaml
---
description: Code style guidelines
alwaysApply: true
---

## Code Style

- Use 2-space indentation
- Prefer const over let
- Use TypeScript strict mode
```

### After → Add to CLAUDE.md

```markdown
# Project Memory

@import .cursor/rules/code-style.mdc

## Additional Guidelines
...
```

**Key changes:**
- No conversion to skill (always-apply rules stay as rules)
- Added `@import` to CLAUDE.md for consistent application

---

## Example 4: Glob-Scoped Rule (No Conversion)

### Before (`.cursor/rules/react-components.mdc`)

```yaml
---
description: React component guidelines
globs: "src/components/**/*.tsx"
alwaysApply: false
---

## React Components

- Use functional components with hooks
- Prop types with TypeScript interfaces
- Co-locate styles with components
```

### After → Keep as Rule

**No conversion needed.** Skills don't support glob-based activation. Keep the rule in place.

---

## Example 5: Background Knowledge Rule

### Before (`.claude/rules/security.md`)

```yaml
---
description: Security best practices
---

## Security Guidelines

1. Always sanitize user input
2. Use parameterized queries
3. Validate file uploads
4. Implement rate limiting
```

### After (`.claude/skills/security-guidelines/SKILL.md`)

```yaml
---
name: security-guidelines
description: |
  Security best practices for the codebase.
  Auto-activates when handling user input, forms, file uploads, 
  database queries, or authentication logic.
---

# Security Guidelines

1. Always sanitize user input
2. Use parameterized queries  
3. Validate file uploads
4. Implement rate limiting

## When to Apply

Consider these guidelines when:
- Handling form submissions
- Writing database queries
- Processing file uploads
- Implementing authentication
- Building API endpoints
```

**Key changes:**
- Created skill directory structure
- Enhanced description with auto-trigger keywords
- Added "When to Apply" section for clarity
- No `disable-model-invocation` (allows auto-activation)

---

## Example 6: Subagent with Hooks

### Before (`.claude/agents/db-validator.md`)

```yaml
---
name: db-validator
description: Validate database queries are read-only
model: haiku
allowed-tools:
  - Bash
  - Read
hooks:
  PreToolUse:
    - matcher: "Bash"
      hooks:
        - type: command
          command: "./scripts/validate-readonly.sh"
---

Execute read-only database queries safely.
Only SELECT statements are permitted.
```

### After (`.claude/skills/db-validator/SKILL.md`)

```yaml
---
name: db-validator
description: |
  Execute and validate database queries ensuring read-only access.
  Use for safe database exploration and data retrieval.
context: fork
model: haiku
allowed-tools: Bash Read
hooks:
  PreToolUse:
    - matcher: "Bash"
      hooks:
        - type: command
          command: "./scripts/validate-readonly.sh"
---

# Database Query Validator

Execute read-only database queries safely.
Only SELECT statements are permitted.

## Safety Rules

1. Only SELECT queries allowed
2. No INSERT, UPDATE, DELETE, DROP
3. Queries validated before execution
4. Results summarized for context efficiency
```

**Key changes:**
- Added `context: fork`
- Simplified model name
- Converted `allowed-tools` to space-delimited format
- Hooks preserved (requires Claude Code 2.1+)

---

## Example 7: Codex Custom Prompt

### Before (`.codex/prompts/format-code.md`)

```yaml
---
name: format-code
description: Format code according to project standards
model: gpt-4o
tools:
  - read
  - write
  - bash
---

Format the provided code files using the project's configured formatters.

1. Detect file type
2. Run appropriate formatter (prettier, black, gofmt, etc.)
3. Report changes made
```

### After (`.claude/skills/format-code/SKILL.md`)

```yaml
---
name: format-code
description: |
  Format code according to project standards.
  Invoke with /format-code <files> to format specific files.
disable-model-invocation: true
argument-hint: <files>
metadata:
  original-model: gpt-4o
allowed-tools: Read Write Bash
---

# Format Code

Format the provided code files using the project's configured formatters.

## Process

1. Detect file type
2. Run appropriate formatter (prettier, black, gofmt, etc.)
3. Report changes made

## Arguments

- `$ARGUMENTS` - Files to format
- `$1`, `$2`, etc. - Individual file paths
```

**Key changes:**
- Standardized tool names (lowercase → PascalCase)
- Preserved original model in metadata
- Added `disable-model-invocation: true` for user invocation
- Added argument documentation

---

## Conversion Summary Template

Use this format when presenting conversion results:

```markdown
## Conversion Summary

### Subagents → Skills (X total)
| Original | New Location | Key Changes |
|----------|--------------|-------------|
| .claude/agents/reviewer.md | .claude/skills/reviewer/SKILL.md | +context:fork, model simplified |

### Commands → Skills (X total)
| Original | New Location | Key Changes |
|----------|--------------|-------------|
| .claude/commands/test.md | .claude/skills/test/SKILL.md | +disable-model-invocation |

### Rules → Various (X total)
| Original | Destination | Action |
|----------|-------------|--------|
| .cursor/rules/style.mdc | CLAUDE.md | @import added |
| .claude/rules/api.md | *(kept)* | Has globs |
| .claude/rules/security.md | .claude/skills/security/SKILL.md | Converted |

### No Action Needed
- .claude/skills/existing/ ✓ (already a skill)
- .cursor/rules/react.mdc ✓ (glob-scoped rule)
```
