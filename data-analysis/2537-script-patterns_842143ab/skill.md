# Script & CLI Patterns Reference

Intelligence for recognizing when a workflow step should be a script, how to design
it as a proper CLI, and how to wire it into a skill. Load before auditing Dimension 4
or scanning for script opportunities during generation.

---

## Signal Patterns: When a Step Should Be a Script

A workflow step is a CLI candidate when any of the following are true. The more signals
present, the stronger the case for scripting.

### Signal 1 — Repeated Generation

The step produces the same structure with different parameters across invocations.
Examples: scaffolding a directory tree, generating a frontmatter block, creating a
boilerplate file from a template. If Claude is re-generating the same code on every
invocation, a parameterized script produces it once and runs reliably thereafter.

*Test:* Would two different users invoking this skill with different inputs cause Claude
to write nearly identical code blocks with just the variable parts swapped? → Script it.

### Signal 2 — Unclear Tool Choice

The step needs to *do something* but no standard Claude tool (Read, Grep, Bash, Edit,
etc.) covers it cleanly without combining multiple tools in a fragile sequence. Example:
"validate frontmatter YAML and report structured errors" requires reading a file, parsing
YAML, and applying rules — awkward as a tool sequence, natural as a script.

*Test:* Does the skill body describe a multi-step procedure that would be done the same
way every time, using tools as primitives? → The procedure is a script waiting to be named.

### Signal 3 — Rigid Input/Output Contract

The step takes a specific input shape (a file path, a name + target directory) and
produces a specific output shape (a scaffolded directory, a JSON report, a validation
result). Rigid contracts are the shape of good CLIs — the interface is clear enough to
parameterize immediately.

*Test:* Can you write the `--help` text for this step right now, without ambiguity?
If yes, it's a CLI. If the args feel unclear, it's still agentic reasoning.

### Signal 4 — Dual-Use Potential

The step would be useful to run independently, outside the skill workflow. Example: a
validation script is useful during skill creation, during repair, and as a standalone
`pre-commit` check. A scaffolding script is useful both when the skill generates a new
artifact and when a user wants to scaffold manually.

*Test:* Would a user want to run this from the terminal directly, without triggering the
full skill? → Design it as a proper CLI from the start, not an internal helper.

### Signal 5 — Consistency Critical

The step must produce identical output for identical inputs — not "similar" output, but
bit-for-bit reproducible results. LLM generation has variance; scripts don't. File
naming conventions, path construction, structural templates — anything where variance
causes downstream breakage should be scripted.

*Test:* Would a subtle difference in output (different field order, different whitespace,
slightly different file name) break something? → Deterministic script, not LLM generation.

---

## CLI Design for Skill Context

A script in a skill directory is also a CLI. Design it to be invoked both by Claude
during a workflow *and* by users from the terminal.

### Interface Design

**Positional arguments** — use for required, ordered inputs where the meaning is
unambiguous from context. Best for 1–2 inputs: `init_skill.py <name> <target-dir>`.

**Named flags** — use for optional inputs, boolean toggles, and anything where the
label clarifies meaning: `--model sonnet`, `--dry-run`, `--output json`.

**Flag for output format** — always add `--output [text|json]` when the script produces
structured data. Claude parses JSON efficiently; humans prefer text. Defaulting to
`text` with `--output json` as the machine-readable mode covers both callers.

**Stdin input** — use when the script is meant to be piped to: `cat file | script.py`.
Useful for transform scripts. Use `sys.stdin.read()` with a flag fallback for file paths.

**Explicit help text** — every script needs `-h`/`--help` output. This is documentation
that Claude reads when deciding how to invoke the script, and that users see when
running it manually. Include: what the script does, each argument/flag with type and
default, and an example invocation.

### Output Conventions

**Stdout for result data** — the primary output goes to stdout. Claude captures stdout.

**Stderr for diagnostic messages** — progress notes, warnings, verbose logging go to
stderr. Claude ignores stderr by default; it doesn't pollute the captured result.

**Exit codes** — `0` for success, `1` for usage errors (wrong args), `2` for runtime
errors (file not found, parse failure). Claude checks exit codes implicitly; a non-zero
exit signals failure and stops the workflow.

**Structured output for multi-field results** — if the script returns more than one
piece of data, output JSON on stdout. A script that returns `{"valid": true, "errors": []}`
is easier for Claude to parse than "Validation passed with 0 errors."

### Script Anatomy (Python template)

```python
#!/usr/bin/env python3
"""
One-line description of what this script does.

Usage:
  script.py <required-arg> [--flag value]

Examples:
  script.py input.yaml --output json
"""

import argparse
import json
import sys

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", help="Description of required input")
    parser.add_argument("--output", choices=["text", "json"], default="text",
                        help="Output format (default: text)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would happen without making changes")
    args = parser.parse_args()

    # Core logic here
    result = process(args.input, dry_run=args.dry_run)

    if args.output == "json":
        print(json.dumps(result))
    else:
        print(format_text(result))

    sys.exit(0 if result["success"] else 1)

if __name__ == "__main__":
    main()
```

---

## Common Script Archetypes

These archetypes cover most script candidates that appear in skill workflows.

### Init — Scaffold a structure

Creates a directory tree or file set from a template. Takes a name and target path;
produces the scaffolded output. Should be idempotent with a `--force` flag for
overwriting, or fail fast on collision by default.

Canonical args: `init.py <name> [target-dir] [--force] [--output json]`

### Validate — Check preconditions

Reads an artifact (file, directory, config), applies a rule set, and reports violations.
Output should be structured (list of `{field, message, severity}` objects). Exit 0 on
clean, exit 1 on violations. Never modifies anything.

Canonical args: `validate.py <path> [--strict] [--output json]`

### Transform — Convert input to output

Takes structured input, applies a deterministic transformation, produces structured
output. The purest CLI form: one input, one output, no side effects unless `--write` is
passed. Use stdin/stdout for pipeline composability.

Canonical args: `transform.py <input-path> [--output-path path] [--dry-run]`

### Package — Assemble an artifact

Collects files or content from multiple sources and assembles a distributable artifact
(zip, tarball, manifest). Should validate inputs before assembling and report what was
included. Dry-run support is valuable here.

Canonical args: `package.py <source-dir> [output-dir] [--dry-run] [--output json]`

### Query — Read state, return structured result

Reads from a data source (DB, file, API) and returns structured data. Never writes.
The primary consumer is Claude reading the output during a skill workflow, but users
should be able to run it for inspection.

Canonical args: `query.py [--filter key=value] [--limit N] [--output json]`

---

## Wiring Scripts into a Skill

A script that isn't referenced in SKILL.md is invisible to Claude.

**In SKILL.md body**, reference each script with:
1. When to invoke it (the trigger condition — which phase, what signals)
2. The exact invocation with relevant flags
3. How to interpret the output (what to do with exit codes, what fields matter)

Example reference pattern:
```markdown
**Validate before proceeding:**
```bash
~/.claude/skills/skill-name/scripts/validate.py "$PATH" --output json
```
Exit 1 = validation failed; parse the `errors` array and report to user before
continuing. Exit 0 = proceed to Phase 3.
```

**Avoid vague references** like "run the validation script if needed" — Claude won't
know which script or when "if needed" applies. State the trigger condition explicitly.

---

## Delegation Pattern

When a skill workflow step is identified as a script candidate, delegate interface
design to the `create-cli` skill rather than designing it ad-hoc. `create-cli` covers
argument structure, help text, output formats, error messages, exit codes, and
config/env precedence in depth.

Invocation pattern from within a skill workflow:
```
Skill: create-cli
Args: "<description of what the script does and what inputs it takes>"
```

The generated CLI spec can then be scaffolded into `scripts/` and referenced from
SKILL.md. This ensures the script is designed for both Claude invocation and direct
user use from the start.
