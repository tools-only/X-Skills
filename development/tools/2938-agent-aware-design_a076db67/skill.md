# Agent-Aware CLI Design — Spec

## Overview

CLI tools built in this workflow are invoked by both AI agents and humans, with agents as the
primary runtime caller. The goal is not to redesign CLI calling conventions — agents are
trained on `--help`, standard flags, and POSIX norms and will use them — but to make every
surface of the existing contract work better when the caller is an LLM. The three pillars are:
token-efficient output, fewer tool calls through good UX, and structured errors that enable
programmatic recovery.

## Key Themes

### Agent-Aware, Not Agent-First

Agents don't need a new interface — they're already fluent with `--help`, flag tables, and
exit codes. What they need is for those surfaces to be precise and low-noise. The skill should
not teach a different CLI paradigm; it should teach authors to be rigorous about the conventions
they already know, and to treat the agent as a peer consumer alongside the human.

### TTY Auto-Detection as the Unifying Mechanism

The output default should be: pretty/human when stdout is a TTY, structured JSON when
piped/non-TTY (which agents always are). This serves both audiences without any required flags.
`--json` and `--human` remain available as explicit overrides. This replaces the current
convention of "add `--json` as an option" — instead, JSON is the automatic default for any
non-interactive caller, including agents.

### Reduce Tool Calls Through CLI Design

Three patterns reduce the number of Bash invocations an agent needs:

- **Compound output:** Operations return enough data to avoid a follow-up call. `create`
  returns the created resource's ID and key fields on stdout. `delete` echoes what was deleted.
  `list` returns full objects, not just IDs, in JSON mode.

- **Rich defaults:** In non-TTY mode, include enough context in a single response that the
  agent rarely needs to call again. Avoid outputs that are "half the answer."

- **Predictable behavior:** Idempotent commands, clear preconditions, consistent flag behavior
  across subcommands. If the agent can predict the outcome, it may skip a verification call
  entirely.

### Consistent Surface for Fast Discovery

Agents scan a CLI's surface area from `--help` output and flag tables. Consistency across
subcommands is what makes that scan cheap:

- Verb-noun subcommand naming, applied uniformly (`create`, `list`, `delete`, `get`).
- Same flag names for the same concepts across all subcommands (`--id`, `--json`, `--force`).
- Same output shape for similar operations (all create-type commands return the same fields).

A compact top-level `--help` that fits in a single screen is the agent's map. It should be
tight enough that the agent reads it once and knows the full surface area.

### Structured Errors with Executable Hints

Error output in non-TTY mode should follow a standard JSON schema on stderr:

```json
{"error": "not_found", "message": "Snapshot 'abc123' does not exist.", "hint": "snapr list --json"}
```

The `hint` field should be an exact command the agent can execute, not prose. This eliminates
a reasoning step — the agent doesn't need to infer what to do, it can run the hint directly.
The three fields are: `error` (machine-readable code, snake_case), `message` (human-readable
sentence), `hint` (optional; exact CLI invocation or `null`).

### Agent Piping Principles (brief)

- Stdout should be stable across versions in non-TTY mode (agents pipe into other tools).
- When a command produces a list, emit one JSON object per line (NDJSON) rather than a
  JSON array — this enables streaming and `jq` piping without buffering the full output.
- Support `-` as stdin/stdout for single-resource commands where it makes sense.
- Never emit ANSI codes, progress spinners, or interactive prompts when stdout is not a TTY.

### Spec Lives in a Skill

The CLI spec produced by this skill typically lives inside a skill body or as an embedded
reference block, not as a standalone file in the user's repo. This means the spec must be
compact enough to fit in an agent's context budget. Redundant sections should be omitted
(not just marked optional), and examples should be dense — demonstrate multiple patterns
in a single invocation rather than one-pattern-per-line.

## Decisions & Positions

- TTY auto-detection is the output mechanism. `--json` is an explicit override, not the
  primary agent affordance.
- Phase 3 of the skill should be rewritten from scratch as unified "CLI Conventions" — not
  "human norms + agent add-on."
- `cli-guidelines.md` should be extended with an "Agent Ergonomics" section (our opinionated
  fork of clig.dev, not a replacement of it).
- The snapr example should be updated to demonstrate all agent-aware patterns: TTY detection,
  structured errors with executable hints, compound output, NDJSON for list commands.
- Flag design stays natural language — no formal schemas. Consistency and precision in the
  spec is sufficient; agents infer well from well-written docs.

## Open Questions

- Should NDJSON (one object per line) be the recommended list format, or should it remain
  a suggestion? NDJSON is strictly more useful for agents but breaks `JSON.parse()` without
  a line-by-line reader.
- Should the error schema require `hint` to be a fully executable command (strict), or allow
  it to be any actionable string (flexible)? Strict is more useful for agents; flexible is
  easier for authors to comply with.

## Constraints & Boundaries

- This is not about new CLI calling conventions. Agents use `--help`, flags, and exit codes
  as-is. The improvements are in output format, error structure, and design discipline.
- Not a formal schema spec (no OpenAPI/JSON Schema for CLI interfaces). Principles and
  strong conventions, not machine-readable contracts.
- Not a token-minimization guide. The goal is fewer tool calls, not shortest possible output.
  Rich compound output often uses *more* tokens per call but fewer calls total.
