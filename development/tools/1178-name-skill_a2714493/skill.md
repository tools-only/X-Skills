---
name: create-cli
description: >
  This skill should be used when the user asks to "design a CLI", "help me design
  command-line flags", "what flags should my tool have", "create a CLI spec",
  "refactor my CLI interface", "design a CLI my agent can call", or wants to design
  command-line UX (args/flags/subcommands/help/output/errors/config) before
  implementation or audit an existing CLI surface for consistency and composability.
argument-hint: "[tool-name and one-line description]"
allowed-tools:
  - AskUserQuestion
  - Read
---

# Create CLI

Design CLI surface area (syntax + behavior), agent-aware, human-friendly.

## Phase 1 — Prepare

Read `${CLAUDE_PLUGIN_ROOT}/skills/create-cli/references/cli-guidelines.md` and
`${CLAUDE_PLUGIN_ROOT}/skills/create-cli/references/language-selection.md`.
Apply cli-guidelines.md as the default CLI rubric; use language-selection.md during Phase 2
to inform the language recommendation.

Proceed when both files are loaded.

## Phase 2 — Clarify

Ask, then proceed with best-guess defaults if user is unsure:

- Command name + one-sentence purpose.
- Primary consumer: agent/LLM, human at a terminal, scripted automation, or mixed.
- Input sources: args vs stdin; files vs URLs; secrets (never via flags).
- Output contract: human text, `--json`, `--plain`, exit codes.
- Interactivity: prompts allowed? need `--no-input`? confirmations for destructive ops?
- Config model: flags/env/config-file; precedence; XDG vs repo-local.
- Language & distribution: ask for the user's preferred implementation language, or offer to
  recommend one. Ask whether a single binary (no runtime needed on target machine) is required,
  or whether a runtime dependency is acceptable. Apply language-selection.md to recommend if
  the user is unsure. Platform: macOS/Linux/Windows.

If an existing CLI spec or tool description is provided, read it first — skip questions already answered by it.

Proceed when answers are confirmed or user is unsure — use best-guess defaults.

## Phase 3 — Conventions

Apply these unless the user says otherwise:

### Output
- Detect TTY: pretty/colored output when stdout is a TTY; structured JSON when piped/non-TTY (agents are always non-TTY). `--json` and `--human` available as explicit overrides.
- List commands: NDJSON (one JSON object per line) in non-TTY mode, not a JSON array — enables streaming and `jq` piping without buffering.
- Primary data to stdout; diagnostics/errors to stderr.
- Suppress ANSI codes, progress spinners, and decorative output in non-TTY mode.

### Errors
- Non-TTY error object on stderr: `{"error": "<snake_case_code>", "message": "...", "hint": "<exact CLI invocation or null>"}` — so agent callers can route recovery logic without parsing free-text stderr. The `hint` field must be an executable command, not prose.
- Exit codes: `0` success, `1` runtime error, `2` invalid usage; add command-specific codes only when genuinely useful.

### Flags
- `-h/--help` always shows help; ignores other args.
- `--version` prints version to stdout.
- Consistent flag names across all subcommands for the same concept (`--id`, `--force`, `--json`) — agents learn the naming pattern once and apply it everywhere without guessing.
- Prompts only when stdin is a TTY; `--no-input` disables prompts.
- Destructive operations: interactive confirmation; non-interactive requires `--force`.
- Respect `NO_COLOR`, `TERM=dumb`; provide `--no-color`.
- Handle Ctrl-C: exit fast; bounded cleanup; crash-only when possible.

### Reduce tool calls
- Compound output: operations return enough data to avoid a follow-up call. `create` returns the new resource's ID and key fields. `delete` echoes what was removed.
- Rich non-TTY defaults: in JSON mode, return full objects not just IDs.
- Idempotent by default: where possible, commands are safe to repeat; document preconditions explicitly — agents rely on safe retries for error recovery without human intervention.

For deeper context on the reasoning behind these conventions, read `${CLAUDE_PLUGIN_ROOT}/skills/create-cli/references/agent-aware-design.md`.

Apply all conventions, then proceed to Phase 4.

## Phase 4 — Deliver

For audits of existing CLIs, produce a gap report (violations + recommended changes) rather than a full spec. For new designs, produce a compact spec the user can implement. Include all relevant sections:

- Command tree + USAGE synopsis.
- Args/flags table (types, defaults, required/optional, examples).
- Subcommand semantics (what each does; idempotence; state changes).
- Output rules: stdout vs stderr; TTY detection; `--json`/`--plain`; `--quiet`/`--verbose`.
- Error + exit code map (top failure modes).
- Safety rules: `--dry-run`, confirmations, `--force`, `--no-input`.
- Config/env rules + precedence (flags > env > project config > user config > system).
- Shell completion story (if relevant): install/discoverability; generation command or bundled scripts.
- 5–10 example invocations (common flows; include piped/stdin examples).

Use this skeleton, dropping irrelevant sections:

0. **Language & distribution**: `Go` · `cobra` · single binary · `goreleaser` for CI
   *(Omit if language was not determined.)*
1. **Name**: `mycmd`
2. **One-liner**: `...`
3. **USAGE**:
   - `mycmd [global flags] <subcommand> [args]`
4. **Subcommands**:
   - `mycmd init ...`
   - `mycmd run ...`
5. **Global flags**:
   - `-h, --help`
   - `--version`
   - `-q, --quiet` / `-v, --verbose` (define exactly)
   - `--human` (force human output when piped) / `--json` (force JSON when TTY) — if TTY auto-detection applies
6. **I/O contract**:
   - stdout:
   - stderr:
7. **Exit codes**:
   - `0` success
   - `1` generic failure
   - `2` invalid usage (parse/validation)
   - (add command-specific codes only when actually useful)
8. **Env/config**:
   - env vars:
   - config file path + precedence:
9. **Examples**:
   - …

See `${CLAUDE_PLUGIN_ROOT}/skills/create-cli/examples/example-cli-spec.md` for a complete worked example.

If the spec is destined for a skill body or CLAUDE.md, omit unused sections entirely (do not mark them "N/A") and limit examples to ≤5 invocations that each demonstrate multiple patterns.

Skill is complete when the spec or gap report is delivered.

## Notes

- Once language is selected (Phase 2), include the idiomatic parsing library in the spec (see language-selection.md). If language remains undetermined, omit the library.
- If the request is "design parameters", do not drift into implementation.
