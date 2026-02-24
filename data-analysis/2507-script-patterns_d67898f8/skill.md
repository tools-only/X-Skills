# Script & CLI Patterns Reference (OpenClaw)

Intelligence for recognizing when a workflow step should be a script, how to design
it as a proper CLI, and how to wire it into an OpenClaw skill. Load before auditing
script opportunities during generation (Phase 2, Step 4).

OpenClaw-specific differences from Claude Code:
- Scripts are invoked via the `exec` tool (not `Bash`)
- Skill-relative paths use `{baseDir}` (not `$CLAUDE_PLUGIN_ROOT`)
- `Glob` and `Grep` tools do not exist — use `exec` + `find`/`rg`/`ls` for filesystem ops
- Tool scoping (e.g., `Bash(git:*)`) is a Claude Code feature — tool policy is gateway config in OpenClaw

---

## Signal Patterns: When a Step Should Be a Script

A workflow step is a CLI candidate when any of the following are true. The more signals
present, the stronger the case for scripting.

### Signal 1 — Repeated Generation

The step produces the same structure with different parameters across invocations.
Examples: scaffolding a directory tree, generating a frontmatter block, creating a
boilerplate file from a template. If the model is re-generating the same code on every
invocation, a parameterized script produces it once and runs reliably thereafter.

*Test:* Would two different users invoking this skill with different inputs cause the model
to write nearly identical code blocks with just the variable parts swapped? → Script it.

### Signal 2 — Unclear Tool Choice

The step needs to *do something* but no standard OpenClaw tool (`read`, `exec`, `write`,
`edit`, etc.) covers it cleanly without combining multiple tools in a fragile sequence.
Example: "validate frontmatter YAML and report structured errors" requires reading a file,
parsing YAML, and applying rules — awkward as a tool sequence, natural as a script.

*Test:* Does the skill body describe a multi-step procedure that would be done the same
way every time, using tools as primitives? → The procedure is a script waiting to be named.

### Signal 3 — Rigid Input/Output Contract

The step takes a specific input shape and produces a specific output shape. Rigid contracts
are the shape of good CLIs — the interface is clear enough to parameterize immediately.

*Test:* Can you write the `--help` text for this step right now, without ambiguity?
If yes, it's a CLI. If the args feel unclear, it's still agentic reasoning.

### Signal 4 — Dual-Use Potential

The step would be useful to run independently, outside the skill workflow. Example: a
validation script is useful during skill creation, during repair, and as a standalone
pre-commit check. A scaffolding script is useful both when the skill generates a new
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

A script in a skill directory is also a CLI. Design it to be invoked both by the model
during a workflow *and* by users from the terminal.

### Interface Design

**Positional arguments** — use for required, ordered inputs where meaning is unambiguous
from context. Best for 1–2 inputs: `init_claw_skill.py <name> <target-dir>`.

**Named flags** — use for optional inputs, boolean toggles, and anything where the label
clarifies meaning: `--dry-run`, `--output json`.

**Flag for output format** — always add `--output [text|json]` when the script produces
structured data. The model parses JSON efficiently; humans prefer text.

**Stdin input** — use when the script is meant to be piped to. Use `sys.stdin.read()` with
a flag fallback for file paths.

**Explicit help text** — every script needs `-h`/`--help`. This is documentation the model
reads when deciding how to invoke the script, and that users see when running it manually.

### Output Conventions

**Stdout for result data** — primary output goes to stdout. The model captures stdout.

**Stderr for diagnostic messages** — progress notes, warnings, verbose logging go to stderr.

**Exit codes** — `0` for success, `1` for usage/validation errors, `2` for runtime errors
(file not found, parse failure).

**Structured output for multi-field results** — if the script returns more than one piece
of data, output JSON on stdout. `{"valid": true, "errors": []}` is easier for the model to
parse than "Validation passed with 0 errors."

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

### Init — Scaffold a structure

Creates a directory tree or file set from a template. Takes a name and target path;
produces the scaffolded output. Should fail fast on collision by default.

Canonical args: `init.py <name> [target-dir] [--output json]`

### Validate — Check preconditions

Reads an artifact (file, directory, config), applies a rule set, and reports violations.
Output should be structured (list of `{field, message, severity}` objects). Exit 0 on
clean, exit 1 on violations. Never modifies anything.

Canonical args: `validate.py <path> [--strict] [--output json]`

### Transform — Convert input to output

Takes structured input, applies a deterministic transformation, produces structured
output. One input, one output, no side effects unless `--write` is passed.

Canonical args: `transform.py <input-path> [--output-path path] [--dry-run]`

### Package — Assemble an artifact

Collects files or content from multiple sources and assembles a distributable artifact.
Should validate inputs before assembling and report what was included.

Canonical args: `package.py <source-dir> [output-dir] [--dry-run] [--output json]`

### Query — Read state, return structured result

Reads from a data source (DB, file, API) and returns structured data. Never writes.

Canonical args: `query.py [--filter key=value] [--limit N] [--output json]`

---

## Wiring Scripts into an OpenClaw Skill

A script that isn't referenced in SKILL.md is invisible to the model.

**In SKILL.md body**, reference each script with:
1. When to invoke it (the trigger condition — which phase, what signals)
2. The exact invocation via `exec` tool, with relevant flags and `{baseDir}` path
3. How to interpret the output (exit codes, which output fields matter)

Example reference pattern:
```markdown
**Validate before proceeding:**

```bash
python3 {baseDir}/scripts/validate_claw_skill.py "$SKILL_DIR" --output json
```

Exit 1 = parse the `errors` array; resolve all `critical` and `major` items before
continuing. Exit 0 = proceed to Phase 3.
```

**Avoid vague references** like "run the validation script if needed" — the model won't
know which script or when "if needed" applies. State the trigger condition explicitly.

Note: Use `exec` tool to run these scripts, not `Bash` (which is not available in OpenClaw
unless the gateway explicitly enables it). Reference paths with `{baseDir}`, not
`$CLAUDE_PLUGIN_ROOT` or hardcoded absolute paths.
