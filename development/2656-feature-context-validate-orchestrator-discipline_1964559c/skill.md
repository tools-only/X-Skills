# Feature Context: Validate and Fix orchestrator-discipline Plugin

## Document Metadata

- **Generated**: 2026-02-20
- **Input Type**: existing_document (groomed context manifest from orchestrator)
- **Source**: Feature request with pre-groomed codebase observations
- **Status**: DISCOVERY_COMPLETE

---

## Original Request

Validate and fix the `orchestrator-discipline` plugin. Grooming has already been done. Identify ambiguities, document use scenarios, and surface questions for resolution.

---

## Feature Summary

The `orchestrator-discipline` plugin enforces delegation discipline in multi-agent Claude Code workflows via two PreToolUse hooks and a rules file. The plugin currently fails `claude plugin validate` due to an unrecognized `rules` field in `plugin.json`. Beyond the validation blocker, there are three substantive gaps: the read-warning hook does not fire when `Grep` targets a directory (only fires when `path` has a source file extension), whether `rules/CLAUDE.md` loads into session context is unverified, and the SKILL.md is missing `user-invocable: true`.

---

## Scope

### In Scope

- Fix the `rules` field in `plugins/orchestrator-discipline/.claude-plugin/plugin.json` to pass `claude plugin validate`
- Determine whether `rules/CLAUDE.md` loads automatically or requires a different mechanism, and implement the correct approach
- Extend `hooks/pre-tool-orchestrator-read-warning.cjs` to fire when `Grep` targets a directory path (not just extension-matched file paths)
- Determine whether `skills/orchestrator-discipline/SKILL.md` needs `user-invocable: true` and add it if so
- Verify hooks fire correctly after fixes
- Version bump per convention (patch fix or minor improvement)

### Out of Scope

- Changing the behavioral semantics of either hook (what they warn about, what they allow)
- Adding new hooks or warning patterns
- Changes to `plugins/orchestrator-discipline/skills/orchestrator-discipline/references/investigation-escalation.md` content
- Marketplace.json changes beyond what the auto-sync pre-commit hook produces automatically

---

## Codebase Research

### Similar Patterns Found

#### Pattern 1: hooks.json at plugin root (not in hooks/ subdirectory)

- **Location**: `plugins/hallucination-detector/.claude-plugin/plugin.json:10` and `plugins/hallucination-detector/hooks.json`
- **Relevance**: The `hallucination-detector` plugin uses `"hooks": "./hooks.json"` pointing to a root-level `hooks.json`, identical to orchestrator-discipline's current approach. This confirms the hooks wiring pattern is valid. The difference is that orchestrator-discipline also has a `hooks/` subdirectory containing the actual JS scripts — the `hooks.json` at root is the config, the `hooks/` directory holds the scripts.
- **Reusable**: Pattern is already in use and consistent.

#### Pattern 2: plugin.json without `rules` field

- **Location**: `plugins/agent-orchestration/.claude-plugin/plugin.json` (all fields), `plugins/hallucination-detector/.claude-plugin/plugin.json`
- **Relevance**: No other plugin in this repository uses a `rules` field. The official plugin.json schema (documented in `plugins/plugin-creator/skills/claude-plugins-reference-2026/SKILL.md:32-52`) lists these component path fields: `commands`, `agents`, `skills`, `hooks`, `mcpServers`, `outputStyles`, `lspServers`. `rules` is absent from the schema.
- **Reusable**: The fix is to remove `rules` from `plugin.json` and determine the correct mechanism for loading the rules file.

#### Pattern 3: `commands` pointing to a `skills/` directory

- **Location**: `plugins/agent-orchestration/.claude-plugin/plugin.json:10` — `"commands": ["./skills/agent-orchestration", "./skills/how-to-delegate"]`
- **Relevance**: orchestrator-discipline currently declares `"commands": ["./skills/orchestrator-discipline"]`. This is the same pattern. However the SKILL.md for orchestrator-discipline lacks `user-invocable: true`, which affects whether the skill surfaces as a slash command.
- **Reusable**: Pattern is valid; the `user-invocable` frontmatter question is separate.

#### Pattern 4: `rules/CLAUDE.md` as a file path convention

- **Location**: No other plugin in the repo has a `rules/` directory. The `orchestrator-discipline` plugin is the only one (`Glob("plugins/*/rules")` returned no results).
- **Relevance**: `rules/CLAUDE.md` is not a recognized plugin component path. In Claude Code, session rules are loaded from `.claude/rules/` at the project or user level — not from plugin directories. There is no `rules` field in the plugin.json schema.
- **Reusable**: N/A — the mechanism is unverified; see Q2 in Questions section.

### Existing Infrastructure

- `plugins/plugin-creator/skills/claude-plugins-reference-2026/SKILL.md` — authoritative plugin.json schema reference, consulted directly
- `plugins/plugin-creator/skills/hooks-io-api/SKILL.md` — hook output format reference (not read in this session; for implementer reference)
- `plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs` — 86-line JS hook, fully read
- `plugins/orchestrator-discipline/hooks/pre-tool-diagnostic-command-gate.cjs` — 110-line JS hook, fully read
- `plugins/orchestrator-discipline/.claude-plugin/plugin.json` — manifest, fully read
- `plugins/orchestrator-discipline/hooks.json` — hook wiring config, fully read
- `plugins/orchestrator-discipline/rules/CLAUDE.md` — 143-line rules file, fully read
- `plugins/orchestrator-discipline/skills/orchestrator-discipline/SKILL.md` — 83-line skill, fully read

### Code References

- `plugins/orchestrator-discipline/.claude-plugin/plugin.json:11` — `"rules": ["./rules"]` — the unrecognized field causing validation failure
- `plugins/orchestrator-discipline/.claude-plugin/plugin.json:12` — `"commands": ["./skills/orchestrator-discipline"]` — skills exposed as commands
- `plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs:15-26` — `SOURCE_FILE_EXTENSIONS` regex and `isSourceOrConfigFile()` — the function that gates hook firing; does not handle directory paths
- `plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs:48-58` — `Grep` path extraction and extension check — gap: a bare directory path like `src/` has no extension, so the check returns false and no warning fires
- `plugins/orchestrator-discipline/skills/orchestrator-discipline/SKILL.md:1-3` — frontmatter has only `description:`, no `user-invocable: true`
- `plugins/plugin-creator/skills/claude-plugins-reference-2026/SKILL.md:32-52` — plugin.json schema, no `rules` field listed
- `plugins/plugin-creator/skills/claude-plugins-reference-2026/SKILL.md:148-153` — default directory auto-loading table; no `rules/` row present

---

## Use Scenarios

### Scenario 1: Plugin Install Fails Validation (Current State)

**Actor**: Developer trying to install the orchestrator-discipline plugin
**Trigger**: Runs `claude plugin validate plugins/orchestrator-discipline/`
**Goal**: Confirm plugin is valid before installation
**Expected Outcome**: Validation passes, plugin can be installed
**Current Outcome**: `✘ Found 1 error: root: Unrecognized key: "rules"` — installation blocked

### Scenario 2: Grep Directory Scan Bypasses Hook (Current Gap)

**Actor**: Orchestrator running as Claude in a session with the plugin installed
**Trigger**: Orchestrator calls `Grep(pattern="class.*Service", path="src/")` to scan a source directory
**Goal**: Hook fires to remind orchestrator to delegate rather than self-investigate
**Current Outcome**: Hook does not fire because `"src/"` has no file extension; `isSourceOrConfigFile("src/")` returns false
**After Fix**: Hook fires with the read-warning because the path targets a source directory

### Scenario 3: Read Warning Fires on Direct File Read (Working)

**Actor**: Orchestrator running as Claude
**Trigger**: Orchestrator calls `Read(file_path="src/cli/commands.py")`
**Goal**: Hook reminds orchestrator to check: "Will I edit this file?"
**Outcome**: Hook fires, injects `<orchestrator-read-warning>` context with the falsifiable test

### Scenario 4: Diagnostic Command Gate Fires (Working)

**Actor**: Orchestrator running as Claude
**Trigger**: Orchestrator calls `Bash(command="uv run ruff check .")`
**Goal**: Hook reminds orchestrator to delegate diagnostic output to an Explore agent
**Outcome**: Hook fires, injects `<orchestrator-diagnostic-warning>` context

### Scenario 5: Rules File Loads Into Session (Unverified)

**Actor**: Any Claude session with plugin installed
**Trigger**: Session starts or plugin enables
**Goal**: Rules in `rules/CLAUDE.md` are loaded into session context as behavioral constraints
**Current Status**: Whether this loading occurs is unverified — the `rules` field is not recognized by the validator

### Scenario 6: User Invokes Skill Manually (Gap)

**Actor**: User in a Claude Code session
**Trigger**: User types `/orchestrator-discipline:orchestrator-discipline`
**Goal**: Load skill guidance on delegation discipline
**Current Status**: SKILL.md lacks `user-invocable: true`; slash command registration status unknown

---

## Gap Analysis

### Identified Gaps

| # | Category | Gap Description | Impact |
|---|----------|-----------------|--------|
| 1 | Scope | `rules` field in plugin.json is not a recognized schema key — causes validation failure | Blocks plugin installation via `claude plugin validate` |
| 2 | Behavior | `pre-tool-orchestrator-read-warning.cjs` does not fire when `Grep` uses a directory path (no extension to match) | Hook misses the most common orchestrator investigation pattern: scanning a directory |
| 3 | Integration | Whether `rules/CLAUDE.md` loads into session context is unverified; no other plugin uses this pattern | Rules content may not reach the model at all |
| 4 | Integration | SKILL.md lacks `user-invocable: true`; slash command registration unverified | Users cannot manually invoke the skill to review discipline guidance |
| 5 | Scope | `commands` field points to `./skills/orchestrator-discipline` — unclear if this is intentional or a schema misuse | The skill may be double-registered or not registered correctly |

---

## Questions Requiring Resolution

### Q1: How should `rules/CLAUDE.md` be delivered to the model?

- **Category**: Integration
- **Gap**: The `rules` field is not in the plugin.json schema. No other plugin in this repo has a `rules/` directory. The validator rejects it. But the rules file contains substantive behavioral constraints that the plugin's value depends on.
- **Question**: Should `rules/CLAUDE.md` be included as part of the skill content in `skills/orchestrator-discipline/SKILL.md`, as a reference file under `skills/orchestrator-discipline/references/`, as a separate CLAUDE.md file loaded via another mechanism, or removed because the hooks alone are sufficient?
- **Options**:
  - A) Merge `rules/CLAUDE.md` content into `skills/orchestrator-discipline/SKILL.md` (skill loads when activated)
  - B) Move `rules/CLAUDE.md` to `skills/orchestrator-discipline/references/rules.md` (referenced from SKILL.md)
  - C) Remove `rules/` entirely — the hooks provide structural enforcement and the rules are duplicated in `~/.claude/CLAUDE.md` already
  - D) Keep `rules/CLAUDE.md` but load it via a different mechanism (e.g., SKILL.md includes it via a markdown link or the skill's context loads it)
- **Why It Matters**: If the rules file is not loaded, the plugin's behavioral-constraint layer is silent. If it is merged into the skill, it only loads when the skill is activated — not passively. If removed, the hook-only enforcement may be sufficient given the CLAUDE.md already contains the same rules.
- **Resolution**: _pending_

### Q2: Should `user-invocable: true` be added to the SKILL.md?

- **Category**: Integration
- **Gap**: SKILL.md at `skills/orchestrator-discipline/SKILL.md:1-3` has no `user-invocable` field. The `commands` entry in `plugin.json` points to the skills directory, which may register the skill as a command. Whether the skill currently appears as `/orchestrator-discipline:orchestrator-discipline` is unverified.
- **Question**: Is manual invocation of this skill by the user a desired use case? The plugin is designed to auto-enforce via hooks. If `user-invocable: true` is added, the skill appears in the slash command menu. If not, it still loads automatically via context.
- **Options**:
  - A) Add `user-invocable: true` — explicit slash command for deliberate activation when reviewing discipline
  - B) Leave as-is — skill loads automatically; no manual invocation needed
- **Why It Matters**: Adding `user-invocable` changes the user-facing interface. The `commands` field in plugin.json may already register it; if so, `user-invocable: true` may be redundant or required to prevent double-registration issues.
- **Resolution**: _pending_

### Q3: What is the correct behavior for `Grep` with a directory path in the hook?

- **Category**: Behavior
- **Gap**: `pre-tool-orchestrator-read-warning.cjs` checks `isSourceOrConfigFile(targetPath)` where `targetPath` is `toolInput.path`. When `path` is a directory like `"src/"` or `"plugins/"`, the extension regex returns false and no warning fires. This is the primary failure mode documented in the groomed context.
- **Question**: When the hook receives a `Grep` call where `path` is a directory, should it always fire (any directory), fire only when the directory name suggests source code (e.g., `src/`, `plugins/`, `tests/`), or fire based on the `pattern` parameter instead?
- **Options**:
  - A) Always fire when `Grep.path` is a directory (most conservative — any directory grep by orchestrator is suspect)
  - B) Fire when `Grep.path` is a directory AND the glob parameter or pattern suggests source file targeting
  - C) Fire when `Grep.path` is a directory that does not match an allowed list (e.g., `plan/`, `docs/`)
  - D) Apply the existing `SOURCE_FILE_EXTENSIONS` check to the `Grep.pattern` parameter as a heuristic (if pattern targets `.py`, `.ts` etc., fire even for directory paths)
- **Why It Matters**: Option A could produce noise for legitimate plan-file or documentation searches. Option D avoids that noise but may miss extension-agnostic searches. The right threshold determines how useful vs. annoying the hook is in practice.
- **Resolution**: _pending_

### Q4: Should the `commands` field in plugin.json be removed?

- **Category**: Scope
- **Gap**: The plugin.json declares both `"skills": ["./skills/orchestrator-discipline"]` and `"commands": ["./skills/orchestrator-discipline"]` pointing to the same path. The schema reference notes that `commands` is the "legacy" path field and skills are the recommended approach. Double-registering the same path may cause no harm, or may cause double-registration.
- **Question**: Should `commands` be removed from plugin.json (leaving only `skills`), or is double-registration intentional to ensure the skill appears both as a skill and as a command?
- **Options**:
  - A) Remove `commands` — skills registration alone is sufficient
  - B) Keep both — belt-and-suspenders to ensure the skill is discoverable
- **Why It Matters**: The official docs say "Skills and commands are now unified" — keeping `commands` is either harmless or creates a duplicate entry. Removing it simplifies the manifest and aligns with current convention.
- **Resolution**: _pending_

---

## Goals (Pending Resolution)

_These goals will be finalized after questions are resolved._

1. `claude plugin validate plugins/orchestrator-discipline/` passes with no errors
2. `pre-tool-orchestrator-read-warning.cjs` fires when `Grep` is called with a directory path targeting source code areas
3. `rules/CLAUDE.md` content reaches the model via a mechanism that is confirmed to work (whether that is skill content, reference file, or removal)
4. `skills/orchestrator-discipline/SKILL.md` has correct frontmatter for intended registration behavior (`user-invocable` set appropriately)
5. Plugin version is bumped to reflect the fixes

---

## Risk Areas

### Risk 1: Rules content silently disappears

If `rules/CLAUDE.md` is simply deleted from the plugin without its content being preserved elsewhere, the behavioral constraint layer of the plugin is lost. The hooks fire non-blocking warnings, but the rules provided the "no exemption categories" and "falsifiable test" language that reinforces hook behavior. If the user's `~/.claude/CLAUDE.md` already contains these rules verbatim, this risk is lower — but that overlap is coincidental, not structural.

**Mitigation**: Verify whether `~/.claude/CLAUDE.md` already includes equivalent rules before deleting `rules/CLAUDE.md`. If it does, deletion is safe. If not, the content must be preserved in the skill or a reference file.

### Risk 2: Hook over-fires after directory-path fix

If the fix for Gap 1 (Grep directory paths) makes the hook fire on all `Grep` calls with a directory path, the orchestrator will see the read-warning even when grepping `plan/` or `docs/` — legitimate orchestrator operations. This could degrade signal quality and cause the warning to be ignored.

**Mitigation**: Define the trigger criteria precisely before implementing. Consider an allowlist of non-source directories.

### Risk 3: `commands` double-registration causes unexpected behavior

Declaring both `skills` and `commands` pointing to the same path in plugin.json is undocumented behavior. It may work, create a duplicate, or silently fail. The fix for the `rules` field may surface this as a secondary validation error.

**Mitigation**: Check whether `claude plugin validate` reports the `commands` + `skills` duplicate after `rules` is removed. Remove `commands` if no functional purpose is identified.

### Risk 4: Hook script executable permissions after caching

The plugin.json schema reference notes that scripts must be executable (`chmod +x`). The JS hooks use `node` via the `command` field in `hooks.json`. Hook files at `hooks/pre-tool-orchestrator-read-warning.cjs` and `hooks/pre-tool-diagnostic-command-gate.cjs` start with `#!/usr/bin/env node` but are invoked as `node "${CLAUDE_PLUGIN_ROOT}/hooks/..."` — so executable permission is not required for the hooks themselves. This risk is low but worth confirming.

### Risk 5: Version bump triggers auto-sync manifest side effects

The pre-commit hook `auto_sync_manifests.py` automatically bumps versions when plugin files change. Depending on which files are modified during the fix, the auto-sync may bump the version independently of any manual version bump. The implementer should let auto-sync handle versioning rather than manually editing version fields to avoid conflicts.

---

## Known Constraints

### Plugin.json Schema Constraints

Per `plugins/plugin-creator/skills/claude-plugins-reference-2026/SKILL.md:32-52`:
- Valid component path fields: `commands`, `agents`, `skills`, `hooks`, `mcpServers`, `outputStyles`, `lspServers`
- `rules` is not a valid field — causes `Unrecognized key: "rules"` validation error
- All paths must start with `./`
- `agents` must be an array of individual file paths (not directory strings)

### Hook JS Conventions in This Repo

Per inspection of both hook files and `plugins/hallucination-detector/hooks.json`:
- Hooks are invoked as `node "${CLAUDE_PLUGIN_ROOT}/hooks/script.js"` (not executed directly)
- Hook scripts read JSON from stdin, write JSON to stdout
- Non-blocking output format: `{ hookSpecificOutput: { hookEventName: "PreToolUse", additionalContext: "..." } }`
- Parse errors result in `{}` output and `process.exit(0)` — hooks never block on error
- `CLAUDE_PLUGIN_ROOT` environment variable is available at hook execution time

### Linting Requirements

Per `.claude/CLAUDE.md` (project instructions):
- JavaScript hooks must use shebangs and executable permissions
- Pre-commit hooks run automatically on commit via `prek`
- Conventional commits require a scope: `fix(orchestrator-discipline): ...`
- The `auto_sync_manifests.py` pre-commit hook will auto-bump plugin version and update marketplace.json

### Skill Name Field Bug (Plugin Skills)

Per `plugins/plugin-creator/CLAUDE.md` (skill name field bug section):
- Plugin skills with an explicit `name:` field in SKILL.md frontmatter do NOT appear as slash commands
- The SKILL.md at `skills/orchestrator-discipline/SKILL.md` currently has no `name:` field — this is correct
- `user-invocable: true` must be set for the skill to appear in the slash command menu

---

## Definition of Done

1. `claude plugin validate plugins/orchestrator-discipline/` exits with no errors
2. `pre-tool-orchestrator-read-warning.cjs` fires (injects `additionalContext`) when tested with a `Grep` call using a directory path (e.g., `Grep(pattern="def ", path="src/")`)
3. `pre-tool-orchestrator-read-warning.cjs` continues to fire for `Grep` calls with source file extension paths
4. `pre-tool-orchestrator-read-warning.cjs` does not fire for `Grep` calls targeting `.md` files or `plan/` directories (no regression on legitimate orchestrator operations)
5. `rules/CLAUDE.md` content is either confirmed to load via a verified mechanism, merged into skill/reference files, or removed with documented rationale
6. SKILL.md frontmatter is updated to reflect correct `user-invocable` value per Q2 resolution
7. `plugin.json` contains no unrecognized fields
8. `uv run prek run --files <changed-files>` passes (markdown linting, pre-commit hooks pass)
9. Commit message follows conventional commits format: `fix(orchestrator-discipline): <description>`

---

## Next Steps

After questions Q1-Q4 are resolved:

1. Update "Resolution" fields in Questions section above
2. Finalize Goals section
3. Proceed to RT-ICA completeness assessment
4. Delegate implementation to `@python3-development:python-cli-architect` (for any Python changes) or direct JS edits via `@general-purpose` agent
5. Verify definition of done checklist
