# Linting Discovery Protocol

Language plugins MUST implement a discovery sequence before executing quality gates. The harness does not assume lint/format/typecheck tools are present or configured.

---

## Discovery Sequence

1. **Git hook** — Check for pre-commit or similar hooks that run quality gates
2. **CI config** — Check for GitHub Actions, GitLab CI, etc. that define gate commands
3. **Project config** — Read pyproject.toml, package.json, etc. for tool config
4. **Fallback** — Use manifest defaults; if tool not found, skip that gate with a warning

---

## Order of Precedence

| Source | Use when |
|--------|----------|
| Project-local override | `.claude/quality-gates.md` or equivalent exists |
| Manifest declaration | Default from language manifest |
| Inferred from file types | No manifest; harness infers from extensions |

---

## Gate Skipping

- **typecheck: (none)** — Skip typecheck gate entirely (non-typed languages)
- **Tool not found** — Skip gate, log warning, continue pipeline
- **Config missing** — Use manifest defaults; warn if tool fails at runtime
