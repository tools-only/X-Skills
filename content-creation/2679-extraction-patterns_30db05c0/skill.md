# Extraction Patterns

How to extract AI-usable knowledge atoms from each user-facing documentation type. Apply these patterns in Phase 1 of the `user-docs-to-ai-skill` workflow.

## Table of Contents

1. [Extraction Atom Format](#extraction-atom-format)
2. [How-To Guides](#how-to-guides)
3. [API References](#api-references)
4. [Tutorials](#tutorials)
5. [Conceptual Explanations](#conceptual-explanations)
6. [Examples and Code Samples](#examples-and-code-samples)
7. [What to Preserve Verbatim](#what-to-preserve-verbatim)
8. [What to Abstract](#what-to-abstract)
9. [Anti-Patterns](#anti-patterns)

---

## Extraction Atom Format

Every extracted piece of knowledge is an atom:

```text
ATOM: <one-sentence fact, constraint, parameter, or pattern>
TYPE: command | parameter | constraint | pattern | error | example
SOURCE: <filename:section-heading>
```

Collect all atoms into a flat list before grouping. Do not filter or discard atoms during extraction — filter during grouping.

---

## How-To Guides

How-to guides are procedure documents. They contain steps, conditions, and outcomes.

**Extract:**

- Each numbered step as a separate `TYPE: pattern` atom
- Conditional branches ("if X, then Y") as `TYPE: constraint` atoms
- Prerequisites listed before steps as `TYPE: constraint` atoms
- Warning/caution blocks as `TYPE: constraint` atoms
- The goal stated in the heading as the atom context

**Example — source text:**

```text
## Configure the output directory

Before running, ensure the target directory exists.

1. Set `output.dir` in your config file.
2. If using relative paths, they resolve from the project root.
3. Run `ty check --output-dir ./results` to verify.

> Warning: existing files in the output directory are overwritten without confirmation.
```

**Extracted atoms:**

```text
ATOM: output.dir must be set in config before running
TYPE: parameter
SOURCE: docs/configuration.md:Configure the output directory

ATOM: target directory must exist before running
TYPE: constraint
SOURCE: docs/configuration.md:Configure the output directory

ATOM: relative paths resolve from project root
TYPE: constraint
SOURCE: docs/configuration.md:Configure the output directory

ATOM: `ty check --output-dir ./results` verifies the output directory configuration
TYPE: command
SOURCE: docs/configuration.md:Configure the output directory

ATOM: existing files in the output directory are overwritten without confirmation
TYPE: constraint
SOURCE: docs/configuration.md:Configure the output directory
```

**Skip:**

- Motivational framing ("This guide will help you...")
- Navigation hints ("See the next section for...")
- Version history context unless it changes current behavior

---

## API References

API references are schema documents. They contain parameters, types, defaults, and constraints.

**Extract:**

- Every parameter name, type, and default as a `TYPE: parameter` atom
- Every enum value list — preserve exactly, do not paraphrase
- Every required vs optional distinction as a `TYPE: constraint` atom
- Every error code or error message pattern as a `TYPE: error` atom
- Return type and shape as a `TYPE: pattern` atom

**Preserve verbatim:** parameter names, type names, enum values, error codes, CLI flag syntax. These are identifiers — paraphrasing causes hallucination downstream.

**Example — source text:**

```text
### `check` command

`ty check [OPTIONS] [PATH]`

Options:
  --output-dir PATH     Write diagnostics to this directory. Default: stdout.
  --format [text|json]  Output format. Default: text.
  --strict              Treat warnings as errors. Default: off.
  --python-version VER  Target Python version. Must be 3.8–3.13.
```

**Extracted atoms:**

```text
ATOM: ty check command syntax is `ty check [OPTIONS] [PATH]`
TYPE: command
SOURCE: docs/cli.md:check command

ATOM: --output-dir PATH writes diagnostics to a directory; default is stdout
TYPE: parameter
SOURCE: docs/cli.md:check command

ATOM: --format accepts [text|json]; default is text
TYPE: parameter
SOURCE: docs/cli.md:check command

ATOM: --strict treats warnings as errors; default is off
TYPE: parameter
SOURCE: docs/cli.md:check command

ATOM: --python-version accepts 3.8 through 3.13 only
TYPE: constraint
SOURCE: docs/cli.md:check command
```

---

## Tutorials

Tutorials are learning documents. They introduce concepts through worked examples. They often contain both conceptual explanation AND procedural steps.

**Extract:**

- The concept being taught as a `TYPE: pattern` atom — what it IS, not why it's interesting
- Each code example as a `TYPE: example` atom — preserve the code verbatim
- Each "why this works" explanation distilled to a single constraint or pattern atom
- Prerequisites listed or implied at the start as `TYPE: constraint` atoms

**Skip:**

- Narrative scaffolding ("In this tutorial, we will learn...")
- Motivational context ("Understanding X will help you Y")
- Analogies and metaphors — extract only the technical fact the analogy explains

---

## Conceptual Explanations

Conceptual docs explain how something works internally. Extract the facts, not the explanation style.

**Extract:**

- Each named concept as a `TYPE: pattern` atom stating what it does
- Relationships between concepts ("X depends on Y", "X replaces Y") as `TYPE: constraint` atoms
- Limitations stated explicitly as `TYPE: constraint` atoms
- Behavior under edge conditions as `TYPE: constraint` atoms

**Skip:**

- Historical background ("This was introduced in version 2.0 because...")
- Design rationale unless it directly constrains usage
- Comparisons to other tools unless they reveal constraints

---

## Examples and Code Samples

Code examples are the highest-fidelity source. They show exact invocation patterns.

**Extract:**

- Every code block as a `TYPE: example` atom — copy verbatim, do not rewrite
- The comment or heading above the code block as the atom description
- If the example demonstrates a constraint ("notice that X must come before Y"), extract that as a separate `TYPE: constraint` atom

**Preserve verbatim:** All code. Never paraphrase code examples — rewriting introduces subtle errors. If the original code has a bug, copy it exactly and add a note in the atom description.

---

## What to Preserve Verbatim

These categories must be copied character-for-character into reference files:

| Category | Reason |
|----------|--------|
| CLI flag names and syntax | Flags are identifiers — paraphrasing causes wrong invocations |
| Parameter names and types | Type names are formal specifications |
| Enum values | Any deviation causes runtime errors |
| Error messages | Error strings are used for pattern matching |
| Code examples | Rewriting introduces bugs |
| Configuration key names | Keys must match exactly |
| File format syntax | Format specs are normative |

---

## What to Abstract

These categories should be distilled, not copied verbatim:

| Category | Abstraction approach |
|----------|---------------------|
| Narrative introductions | Drop entirely — extract only the technical claim |
| Step-by-step prose | Convert to `ATOM: do X to achieve Y` |
| Analogies | Extract the constraint the analogy illustrates |
| Repetitive examples | Extract one canonical example; note that others exist |
| Version history | Extract only if current behavior differs from prior |
| FAQ sections | Treat each Q/A pair as a constraint atom |

---

## Anti-Patterns

**Summarizing code instead of copying it:**

```text
# WRONG
ATOM: Run the check command with appropriate flags
TYPE: command

# CORRECT
ATOM: `ty check --strict --format json ./src` runs strict checking with JSON output
TYPE: command
SOURCE: docs/quickstart.md:Running checks
```

**Merging multiple constraints into one atom:**

```text
# WRONG
ATOM: The config file must be named ty.toml and placed in the project root and must include the python-version field
TYPE: constraint

# CORRECT — three separate atoms
ATOM: Config file must be named ty.toml
TYPE: constraint
SOURCE: docs/configuration.md:Config file location

ATOM: Config file must be placed in the project root
TYPE: constraint
SOURCE: docs/configuration.md:Config file location

ATOM: Config file must include the python-version field
TYPE: constraint
SOURCE: docs/configuration.md:Required fields
```

**Extracting motivational framing as knowledge:**

```text
# WRONG
ATOM: ty is a fast Python type checker designed for modern workflows
TYPE: pattern

# CORRECT — skip marketing copy; extract the technical fact
ATOM: ty checks Python type annotations without running the code
TYPE: pattern
SOURCE: docs/overview.md:What ty does
```

**Paraphrasing error messages:**

```text
# WRONG
ATOM: ty reports an error when the config file is missing
TYPE: error

# CORRECT — copy the exact message
ATOM: ty emits "error: no ty.toml found in current directory or any parent directory" when no config file exists
TYPE: error
SOURCE: docs/errors.md:Configuration errors
```
