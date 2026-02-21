---
description: "Fix and validate the orchestrator-discipline plugin — blocker, hook coverage, SKILL.md frontmatter, and full validation suite"
version: "1.0"
tasks:
  - T1: Fix plugin.json validation failure (BLOCKER)
  - T2: Fix Grep directory path coverage in hook
  - T3: Fix SKILL.md frontmatter and linting
  - T4: Full validation suite (depends on T1 + T2 + T3)
task_exports:
  enabled: false
  directory: "TASK"
---

# Task Plan: Validate and Fix orchestrator-discipline Plugin

## Context

The `orchestrator-discipline` plugin enforces delegation discipline via two PreToolUse hooks and a
rules file. It currently fails `claude plugin validate` due to an unrecognized `rules` field in
`plugin.json`. Three additional gaps exist: directory-targeted `Grep` bypasses the read-warning
hook, `rules/CLAUDE.md` loading is unverified, and `SKILL.md` is missing `user-invocable: true`.

Source documents:

- `plan/feature-context-validate-orchestrator-discipline.md`
- `plan/codebase/orchestrator-discipline-patterns.md`

Plugin root: `plugins/orchestrator-discipline/`

---

## Context Manifest

### Key Decisions Already Made (Do Not Re-Investigate)

The architecture spec (`plan/architect-validate-orchestrator-discipline.md`) resolved all four open questions from the feature context doc. These decisions are final — do not reopen them.

**Decision 1 — rules/CLAUDE.md delivery mechanism**: Move content verbatim to a new file at `plugins/orchestrator-discipline/CLAUDE.md` (plugin root). This is the same pattern used by `plugins/development-harness/CLAUDE.md` and `plugins/plugin-creator/CLAUDE.md`. The `rules/` directory and `rules/CLAUDE.md` are to be deleted after content is confirmed moved. Option B (merge into SKILL.md body) was rejected because SKILL.md content loads on demand rather than passively. Option C (references file) was also rejected.

**Decision 2 — Grep directory path coverage**: Add `isDirectory()` + `looksLikeDirectory()` to the Grep branch in `pre-tool-orchestrator-read-warning.cjs`. Fire the warning for ALL directory-targeted Grep calls — no allowlist, no pattern heuristic. The hook is non-blocking so the noise cost is acceptable.

**Decision 3 — SKILL.md `user-invocable` field**: Add `user-invocable: true` to SKILL.md frontmatter. Do NOT add a `name:` field. The absence of `name:` is intentional and correct — adding it would suppress slash command registration due to a confirmed Claude Code v2.1.23 bug (documented in `plugins/plugin-creator/CLAUDE.md`, "CRITICAL: Skill Name Field Bug" section).

**Decision 4 — Remove `commands` field from plugin.json**: The `commands` field is the legacy registration path. Both `skills` and `commands` currently point to `./skills/orchestrator-discipline`. Remove `commands`; keep only `skills`. After T1 completes, `plugin.json` should contain exactly: `name`, `skills`, `hooks`.

---

### Current State of Each File to Be Modified

**`plugins/orchestrator-discipline/.claude-plugin/plugin.json`** — current contents (13 lines):

```json
{
  "name": "orchestrator-discipline",
  "description": "Enforces orchestrator context window discipline...",
  "version": "1.3.2",
  "author": { "name": "Jamie Nelson", "url": "https://github.com/bitflight-devops" },
  "skills": ["./skills/orchestrator-discipline"],
  "hooks": "./hooks.json",
  "rules": ["./rules"],
  "commands": ["./skills/orchestrator-discipline"]
}
```

Target state after T1: remove `"rules"` and `"commands"` lines. Do NOT touch `version` — `auto_sync_manifests.py` pre-commit hook handles versioning automatically.

**`plugins/orchestrator-discipline/rules/CLAUDE.md`** — 143-line behavioral constraints file. Content covers: Context Window Read Constraints (permitted/prohibited), Delegation Constraints, Investigation Escalation Anti-Pattern (with mermaid diagram), Agent Output Polling Anti-Pattern (with mermaid diagram), Diagnostic Commands list, Epistemic Identity Scope. Key string that must survive in destination: `"no exemption categories"` (used in T1 and T4 verification).

**`plugins/orchestrator-discipline/CLAUDE.md`** — does NOT exist yet. T1 creates it by copying `rules/CLAUDE.md` content verbatim.

**`plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs`** — 86-line CommonJS Node.js hook. Current structure:

- Lines 15-16: `SOURCE_FILE_EXTENSIONS` regex and `TEST_PATH_PATTERN` regex constants.
- Lines 23-26: `isSourceOrConfigFile(filePath)` function — returns true if path matches extension regex OR test path pattern.
- Lines 28-58: stdin read loop, JSON parse, tool name dispatch. For `Read`, reads `toolInput.file_path`. For `Grep`, reads `toolInput.path`. Falls through to `isSourceOrConfigFile(targetPath)` check — this is the gap (directory paths like `"src/"` return false).
- Lines 60-84: Produces `{ hookSpecificOutput: { hookEventName: "PreToolUse", additionalContext: "<orchestrator-read-warning>...</orchestrator-read-warning>" } }` and exits 0.
- No `require('node:fs')` import exists — must be added in T2.
- Uses CommonJS (`require()`), not ESM (`import`). Must stay CommonJS.

**`plugins/orchestrator-discipline/skills/orchestrator-discipline/SKILL.md`** — 83-line skill. Current frontmatter (lines 1-3):

```yaml
---
description: Orchestrator context window discipline enforcement. Prevents...
---
```

Target state after T3: add `user-invocable: true` on line 3, before closing `---`. No `name:` field. Body content unchanged.

**`plugins/orchestrator-discipline/hooks.json`** — does NOT need modification. It wires two PreToolUse hooks: `Read|Grep` matcher to `pre-tool-orchestrator-read-warning.cjs` and `Bash` matcher to `pre-tool-diagnostic-command-gate.cjs`, both invoked as `node "${CLAUDE_PLUGIN_ROOT}/hooks/..."`.

---

### Hook Architecture — How the Read Warning Hook Works

When Claude Code fires a PreToolUse event for `Read` or `Grep`, it pipes a JSON object to the hook's stdin in the form:

```json
{ "tool_name": "Grep", "tool_input": { "path": "src/", "pattern": "foo" } }
```

The hook reads all stdin, parses the JSON, checks `tool_name`, extracts `targetPath` (from `tool_input.file_path` for Read, from `tool_input.path` for Grep), then calls `isSourceOrConfigFile(targetPath)`. If that returns true, the hook writes a `hookSpecificOutput.additionalContext` JSON block. If false, it writes `{}`. Both paths call `process.exit(0)` — the hook is strictly non-blocking.

The gap: `isSourceOrConfigFile("src/")` returns false because neither the extension regex nor the test-path pattern match a bare directory string. The fix adds two new functions (`isDirectory` and `looksLikeDirectory`) and OR-extends the `shouldWarn` condition in the Grep branch only.

The exact logic change to implement (from architect spec):

```text
// Add at top of file alongside existing requires:
const fs = require('node:fs');

// Add two new helper functions after isSourceOrConfigFile():

function isDirectory(targetPath) {
  try { return fs.statSync(targetPath).isDirectory(); } catch { return false; }
}

function looksLikeDirectory(targetPath) {
  if (!targetPath) return false;
  if (targetPath.endsWith('/') || targetPath.endsWith('\\')) return true;
  const lastSegment = targetPath.split(/[\\/]/).filter(Boolean).pop() || '';
  return !lastSegment.includes('.');
}

// Change Grep branch from:
shouldWarn = isSourceOrConfigFile(targetPath) || isTestPath(targetPath)

// To:
shouldWarn = isSourceOrConfigFile(targetPath)
          || isTestPath(targetPath)
          || isDirectory(targetPath)
          || looksLikeDirectory(targetPath)
```

Note: `BACKLOG.md` has a `.md` extension — it does NOT match `SOURCE_FILE_EXTENSIONS`, is NOT a test path, `fs.statSync("BACKLOG.md").isDirectory()` returns false, and it does NOT pass `looksLikeDirectory` (it has `.` in the last segment). So `BACKLOG.md` will correctly produce `{}` output.

---

### `user-invocable: true` Syntax — Confirmed Pattern

From `plugins/plugin-creator/skills/skill-creator/SKILL.md` lines 1-5 (read directly):

```yaml
---
description: Guide for creating effective skills...
license: Complete terms in LICENSE.txt
user-invocable: true
---
```

The field name is `user-invocable` (hyphenated), value is YAML boolean `true` (not the string `"true"`). Place after `description:`, before the closing `---`.

---

### Plugin.json Schema — What Is and Is Not Allowed

Per `plugins/plugin-creator/skills/claude-plugins-reference-2026/SKILL.md:32-52` and confirmed in `plugins/plugin-creator/CLAUDE.md` "Plugin.json Requirements" section:

Valid component path fields: `commands`, `agents`, `skills`, `hooks`, `mcpServers`, `outputStyles`, `lspServers`.

The field `rules` is NOT in the schema — this is the direct cause of the `Unrecognized key: "rules"` validation error.

The `commands` field is the legacy registration path. Per the CLAUDE.md: "Skills and slash commands are now unified — they are the same system." The canonical registration after the fix is `skills` in plugin.json plus `user-invocable: true` in SKILL.md frontmatter.

---

### Plugin Root CLAUDE.md Pattern — Confirmed

Two existing plugins use root-level `CLAUDE.md` as the mechanism for delivering session-level behavioral context:

- `plugins/development-harness/CLAUDE.md` — confirmed present (directory listing shows it)
- `plugins/plugin-creator/CLAUDE.md` — confirmed present (read in full above)

`plugins/orchestrator-discipline/` currently has NO root-level `CLAUDE.md`. The `rules/CLAUDE.md` at `plugins/orchestrator-discipline/rules/CLAUDE.md` is the content to be moved.

---

### Validation Tools

Two validators are used in T4:

1. `claude plugin validate plugins/orchestrator-discipline/` — official Claude Code CLI validator. Checks: plugin.json schema compliance, `name` field present and kebab-case, all paths start with `./`, referenced files exist. Current failure: `✘ Found 1 error: root: Unrecognized key: "rules"`.

2. `./plugins/plugin-creator/scripts/plugin_validator.py plugins/orchestrator-discipline/` — repo-internal validator. Invoked as `uv run plugins/plugin-creator/scripts/plugin_validator.py plugins/orchestrator-discipline/`. Validates frontmatter, skill complexity (token-based), internal links, cross-references. Use `--verbose` for detail; `--fix` for auto-remediation of known issues (removes `name:` fields from plugin skills automatically).

Linting tool: `uv run prek run --files <path>` — runs all configured pre-commit hooks against specific files. Used to confirm markdown formatting compliance.

---

### Critical Constraints for All Tasks

- Do NOT manually edit `version` fields in `plugin.json` or `.claude-plugin/marketplace.json` — `auto_sync_manifests.py` pre-commit hook handles versioning on commit.
- Do NOT run `git commit` or `git add` — leave staging for human review.
- Do NOT use ESM `import` syntax in the JS hook — file uses CommonJS `require()`.
- Do NOT add external npm dependencies to the hook — use only `node:` stdlib modules.
- Do NOT add a `name:` field to SKILL.md — doing so suppresses slash command registration.
- Do NOT change the `additionalContext` XML content or tag name in the hook.
- Conventional commit scope for this plugin: `fix(orchestrator-discipline): ...`

---

### File Paths for Implementation

| File | Task | Action |
|------|------|--------|
| `plugins/orchestrator-discipline/.claude-plugin/plugin.json` | T1 | Remove `rules` and `commands` fields |
| `plugins/orchestrator-discipline/CLAUDE.md` | T1 | Create — copy verbatim content from `rules/CLAUDE.md` |
| `plugins/orchestrator-discipline/rules/CLAUDE.md` | T1 | Document for deletion (do not delete in T1 — leave for human) |
| `plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs` | T2 | Add `fs` require + two helper functions + extend Grep shouldWarn |
| `plugins/orchestrator-discipline/skills/orchestrator-discipline/SKILL.md` | T3 | Add `user-invocable: true` to frontmatter only |

Reference files (read-only):

- `plugins/plugin-creator/skills/claude-plugins-reference-2026/SKILL.md` — lines 32-52 (schema), lines 148-153 (auto-loading table)
- `plugins/plugin-creator/CLAUDE.md` — skill name field bug, plugin.json requirements
- `plugins/orchestrator-discipline/hooks/pre-tool-diagnostic-command-gate.cjs` — do NOT modify (out of scope)
- `plugins/orchestrator-discipline/skills/orchestrator-discipline/references/investigation-escalation.md` — do NOT modify (out of scope)

---

### What NOT to Re-Investigate

The following questions were resolved by the architect doc and do NOT need re-investigation:

- Whether `rules` is a valid plugin.json field — it is NOT (confirmed schema gap).
- Which delivery mechanism for rules content — plugin root `CLAUDE.md` is chosen (Decision 1).
- Whether to keep `commands` field — remove it (Decision 4).
- Whether `user-invocable: true` is needed — yes (Decision 3).
- Whether `name:` should be added — NO, confirmed bug means omitting `name:` is correct.
- Whether `fs.statSync` is safe in CommonJS Node.js hooks — yes, with try/catch (confirmed pattern).
- Whether `BACKLOG.md` would trigger the new directory detection — confirmed NO (has `.md` extension, `isDirectory()` returns false, `looksLikeDirectory` returns false).

---

## Dependency Graph

```
T1 (plugin.json fix)    ─┐
T2 (Grep dir coverage)  ─┼─→ T4 (Full validation suite)
T3 (SKILL.md frontmatter)┘
```

T1, T2, and T3 are independent and can execute in parallel.
T4 depends on T1, T2, and T3 all completing successfully.

---

## Priority Ordering

| Task | Priority | Complexity | Depends On     | Parallelize With |
|------|----------|------------|----------------|------------------|
| T1   | High     | Low        | None           | T2, T3           |
| T2   | High     | Medium     | None           | T1, T3           |
| T3   | Medium   | Low        | None           | T1, T2           |
| T4   | High     | Medium     | T1 + T2 + T3   | None             |

---

## Tasks

---

### T1: Fix plugin.json Validation Failure

```yaml
---
task: T1
title: Fix plugin.json validation failure (BLOCKER)
status: not-started
agent: general-purpose
dependencies: []
priority: 1
complexity: low
accuracy-risk: medium
parallelize-with: [T2, T3]
reason: T1, T2, T3 touch different files — plugin.json, a JS hook, and SKILL.md respectively. No file conflicts.
handoff: "Report: final plugin.json contents, output of `claude plugin validate`, disposition of rules/CLAUDE.md content (merged/moved/kept), and whether `commands` field was removed."
---
```

#### Context

`plugins/orchestrator-discipline/.claude-plugin/plugin.json` currently contains a `"rules"` key
that is not in the Claude Code plugin.json schema. This causes:

```
✘ Found 1 error: root: Unrecognized key: "rules"
```

The schema (per `plugins/plugin-creator/skills/claude-plugins-reference-2026/SKILL.md:32-52`)
allows only: `commands`, `agents`, `skills`, `hooks`, `mcpServers`, `outputStyles`, `lspServers`.

The `rules` field points to `./rules`, which contains
`plugins/orchestrator-discipline/rules/CLAUDE.md` — a 143-line behavioral constraints file. The
plugin also declares both `"skills"` and `"commands"` pointing to the same
`./skills/orchestrator-discipline` path, which may cause a double-registration issue surfaced after
the `rules` field is removed.

No other plugin in this repository uses a `rules` field. Alternate patterns for surfacing rules
content: plugin root `CLAUDE.md`, inline in `SKILL.md`, or `skills/.../references/` file.

#### Objective

Remove the invalid `rules` field from `plugin.json`, preserve `rules/CLAUDE.md` content via a
verified loading mechanism, and make `claude plugin validate plugins/orchestrator-discipline/` exit
0 with no errors.

#### Required Inputs

- `plugins/orchestrator-discipline/.claude-plugin/plugin.json` — current manifest
- `plugins/plugin-creator/skills/claude-plugins-reference-2026/SKILL.md` — authoritative schema
  reference (lines 32-52 for component path fields, lines 148-153 for auto-loading table)
- `plugins/orchestrator-discipline/rules/CLAUDE.md` — content to be preserved
- `plugins/orchestrator-discipline/skills/orchestrator-discipline/SKILL.md` — may need update if
  content is merged here
- `plugins/plugin-creator/CLAUDE.md` — skill name field bug notes (affects whether `commands`
  field is needed)
- Assumption: `claude` CLI is available on PATH. Confirm with `which claude`.

#### Requirements

1. Remove `"rules": ["./rules"]` from `plugin.json`.
2. Move the content of `./rules/CLAUDE.md` to `plugins/orchestrator-discipline/CLAUDE.md` (plugin
   root). This is the confirmed pattern used by `plugins/development-harness/` and
   `plugins/plugin-creator/` — root-level `CLAUDE.md` is auto-loaded as plugin documentation
   context. Do not choose option B (SKILL.md body) or option C (references/rules.md); the
   architecture spec (Decision 1) has already evaluated and rejected those options.
3. Implement whichever option from requirement 2 is chosen. Do NOT delete `rules/CLAUDE.md`
   content without confirming it exists in the destination.
4. Assess whether `"commands": ["./skills/orchestrator-discipline"]` should be removed. Read
   `plugins/plugin-creator/CLAUDE.md` to confirm whether `commands` + `skills` pointing to the
   same path causes double-registration. Remove `commands` if it is redundant and the schema
   recommends `skills` alone.
5. Run `claude plugin validate plugins/orchestrator-discipline/` and verify exit code is 0.

#### Constraints

- Do NOT change the behavioral semantics of the hooks (what they warn about, what they allow).
- Do NOT modify `hooks/pre-tool-orchestrator-read-warning.cjs` or
  `hooks/pre-tool-diagnostic-command-gate.cjs` — those are T2's scope.
- Do NOT modify `skills/orchestrator-discipline/SKILL.md` frontmatter fields — that is T3's scope.
  You may append content to the body if merging rules, but do not touch frontmatter.
- Do NOT manually bump the version in `plugin.json` — the `auto_sync_manifests.py` pre-commit hook
  handles versioning automatically on commit.
- Do NOT run `git commit` — leave staging for the human to review.

#### Expected Outputs

- `plugins/orchestrator-discipline/.claude-plugin/plugin.json` — modified (no `rules` field, and
  possibly no `commands` field if determined redundant)
- One of:
  - `plugins/orchestrator-discipline/CLAUDE.md` created (if option A chosen), OR
  - `plugins/orchestrator-discipline/skills/orchestrator-discipline/SKILL.md` body updated (if
    option B chosen), OR
  - `plugins/orchestrator-discipline/skills/orchestrator-discipline/references/rules.md` created
    (if option C chosen)
- If `rules/CLAUDE.md` is superseded, document that it can be removed (but do NOT delete it in
  this task — leave deletion decision to human review).

#### Acceptance Criteria

1. `claude plugin validate plugins/orchestrator-discipline/` exits 0 with no errors printed.
2. `plugins/orchestrator-discipline/.claude-plugin/plugin.json` contains no `rules` key when
   inspected with `cat`.
3. The full text of the original `rules/CLAUDE.md` is reachable via one of the approved mechanisms
   (plugin root CLAUDE.md, SKILL.md body, or references/rules.md). Verify by reading the
   destination file and confirming the "no exemption categories" and "falsifiable test" language is
   present.
4. If `commands` field was removed: `plugin.json` contains no `commands` key. If kept: document
   why in handoff.
5. `node -e "JSON.parse(require('fs').readFileSync('plugins/orchestrator-discipline/.claude-plugin/plugin.json','utf8'))"` exits 0 (valid JSON).

#### Verification Steps

1. `claude plugin validate plugins/orchestrator-discipline/` — confirm output shows no errors and
   exit code is 0.
2. `cat plugins/orchestrator-discipline/.claude-plugin/plugin.json` — visually confirm no `rules`
   key present.
3. Read the destination file for `rules/CLAUDE.md` content. Grep for the string
   "no exemption categories" to confirm content survived.
4. `node -e "JSON.parse(require('fs').readFileSync('plugins/orchestrator-discipline/.claude-plugin/plugin.json','utf8'))"` — confirm valid JSON.

#### CoVe Checks

Key claims to verify:

- Claim: `rules` is not in the plugin.json schema.
- Claim: Plugin root `CLAUDE.md` is auto-loaded as plugin documentation context (per how other
  plugins deliver rule content).
- Claim: `commands` and `skills` pointing to the same path may cause double-registration.

Verification questions:

1. Does `plugins/plugin-creator/skills/claude-plugins-reference-2026/SKILL.md` lines 32-52 list
   `rules` as a valid component path field? (Expected: No.)
2. Does `plugins/plugin-creator/skills/claude-plugins-reference-2026/SKILL.md` lines 148-153
   (default directory auto-loading table) include a `rules/` or root `CLAUDE.md` row? (Determines
   which delivery mechanism is verified.)
3. Does `plugins/plugin-creator/CLAUDE.md` document whether `commands` + `skills` both pointing to
   the same directory path causes a schema error or harmless redundancy?

Evidence to collect:

- Read lines 32-52 and 148-153 of
  `plugins/plugin-creator/skills/claude-plugins-reference-2026/SKILL.md`.
- Read `plugins/plugin-creator/CLAUDE.md` for the skill name field bug and `commands` vs `skills`
  guidance.
- Read `plugins/development-harness/CLAUDE.md` or `plugins/agent-orchestration/CLAUDE.md` as an
  example of plugin root `CLAUDE.md` usage (confirm the pattern exists).

Revision rule:

If any CoVe check reveals the schema allows `rules`, or that root `CLAUDE.md` is NOT auto-loaded,
revise the delivery mechanism choice and state what changed in the handoff.

---

### T2: Fix Grep Directory Path Coverage in Hook

```yaml
---
task: T2
title: Fix Grep directory path coverage in pre-tool-orchestrator-read-warning.cjs
status: not-started
agent: general-purpose
dependencies: []
priority: 1
complexity: medium
accuracy-risk: medium
parallelize-with: [T1, T3]
reason: T2 modifies only the JS hook file. T1 touches plugin.json and potentially rules/CLAUDE.md. T3 touches SKILL.md. No overlap.
handoff: "Report: the exact code change made (diff or before/after), output of `node --check`, test evidence for all four Grep scenarios in acceptance criteria, and whether BACKLOG.md Grep calls were determined to trigger or not (with rationale documented)."
---
```

#### Context

`plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs` (86 lines) fires a
`<orchestrator-read-warning>` advisory for `Read` and `Grep` tool calls when the target path has a
source file extension. The current gate function:

```javascript
const SOURCE_FILE_EXTENSIONS =
  /\.(py|toml|yaml|yml|js|ts|jsx|tsx|json|cfg|ini|env|sh|bash|go|rs|rb|java|c|cpp|h|hpp)$/i;

function isSourceOrConfigFile(filePath) {
  return SOURCE_FILE_EXTENSIONS.test(filePath) || TEST_PATH_PATTERN.test(filePath);
}
```

When an orchestrator calls `Grep(pattern="class.*Service", path="src/")`, `targetPath` is `"src/"`.
`isSourceOrConfigFile("src/")` returns false — no extension match, no test path match — and the
hook silently passes. This is the primary failure mode: directory-targeted Grep by an orchestrator
escapes the warning.

The codebase patterns doc (`plan/codebase/orchestrator-discipline-patterns.md`, section 5.2)
provides two detection approaches: `fs.statSync` (filesystem check) and pattern heuristics
(trailing slash or no dot in last segment).

The hook must NOT crash when `path` is absent from `tool_input` (current code already guards with
`toolInput.path || ''`).

#### Objective

Extend the hook so that `Grep` calls where `path` is a directory (detected via heuristic or
`fs.statSync`) also trigger the `<orchestrator-read-warning>` advisory, without introducing
crashes or false-positive noise on legitimate markdown/plan-file searches.

#### Required Inputs

- `plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs` — full file to
  read before editing
- `plan/codebase/orchestrator-discipline-patterns.md` sections 4 and 5 — gap analysis and
  detection options
- `plan/feature-context-validate-orchestrator-discipline.md` Q3 (lines 183-194) — options A–D for
  directory detection scope

#### Requirements

1. Read the full current hook file before making any changes.
2. Add directory detection to the `Grep` path-checking logic. Implement **Option A from Q3** (fire
   for any directory path) combined with a heuristic guard: treat a path as a directory if
   `fs.statSync(path).isDirectory()` returns true, OR if the path ends with `/`, OR if the last
   segment of the path contains no `.` character. Wrap `fs.statSync` in a try/catch that returns
   false on any error (path does not exist or is inaccessible).
3. `BACKLOG.md` is a legitimate orchestrator read target (tracking file, not source code). Grep
   calls with `path: "BACKLOG.md"` must NOT trigger the hook. Confirm that `BACKLOG.md` does not
   match the new directory detection logic (it will not — it has a `.md` extension and
   `isDirectory()` will return false). Document in a code comment that `.md` files are excluded by
   the extension check already.
4. Do not alter the existing behavior for `Read` calls or for `Grep` calls with source file
   extension paths — only extend the `Grep` directory detection branch.
5. Preserve the existing error-safety pattern: if `path` is absent or empty, the hook passes
   through silently (return `{}`).

#### Constraints

- Do NOT change the `additionalContext` message content or XML tag name — only the trigger logic.
- Do NOT change how the hook handles `Read` tool calls.
- Do NOT add any `console.error` or `console.warn` calls — the hook writes only to stdout.
- Do NOT use ES module syntax (`import`) — the file uses CommonJS `require()`.
- Do NOT add external npm dependencies — use only `node:` stdlib modules (`node:fs`, `node:path`).
- The hook must still exit via `process.exit(0)` in all code paths.
- `node --check` must pass on the modified file.

#### Expected Outputs

- `plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs` — modified with
  directory detection logic added to the Grep branch

#### Acceptance Criteria

1. `node --check plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs`
   exits 0 (no syntax errors).
2. A test script that simulates `Grep({ pattern: "foo", path: "src/" })` via stdin produces JSON
   output containing `"additionalContext"`. (Simulate by piping JSON to the script via
   `echo '{"tool_name":"Grep","tool_input":{"path":"src/","pattern":"foo"}}' | node
   plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs`.)
3. A test with `path: "plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs"` (an
   existing `.js` file) produces `"additionalContext"` (existing behavior preserved).
4. A test with `path: "BACKLOG.md"` produces `{}` (no warning — markdown file excluded by
   extension check; not a source file extension, not a directory).
5. A test with `path: ""` (empty string) produces `{}` without crashing.

#### Verification Steps

1. `node --check plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs`
2. `echo '{"tool_name":"Grep","tool_input":{"path":"src/","pattern":"foo"}}' | node plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs` — output must contain `additionalContext`.
3. `echo '{"tool_name":"Grep","tool_input":{"path":"plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs","pattern":"foo"}}' | node plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs` — output must contain `additionalContext`.
4. `echo '{"tool_name":"Grep","tool_input":{"path":"BACKLOG.md","pattern":"foo"}}' | node plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs` — output must be `{}`.
5. `echo '{"tool_name":"Grep","tool_input":{"path":"","pattern":"foo"}}' | node plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs` — output must be `{}` and script must not throw.

#### CoVe Checks

Key claims to verify:

- Claim: `tool_input.path` is the correct field name for the Grep path argument in the hook's
  stdin JSON.
- Claim: `fs.statSync(path).isDirectory()` with a try/catch is safe to call from a Node.js hook
  in the repo's CommonJS + `node:fs` convention.
- Claim: `BACKLOG.md` will not trigger the new directory detection logic.

Verification questions:

1. Does the current hook file read Grep path as `toolInput.path` (not `toolInput.file_path` or
   another field)? Read `plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs`
   lines 48-51 to confirm.
2. Does `BACKLOG.md` end with `/`? Does it have no `.` in the last segment? (Both are false — it
   ends in `.md`.) Confirm that neither heuristic matches it and `fs.statSync("BACKLOG.md").isDirectory()` returns false (it is a file).
3. Does any other hook in this repo use `fs.statSync` from `node:fs`? If not, is there a
   documented reason to avoid it? Check `plugins/*/hooks/*.js` for filesystem access patterns.

Evidence to collect:

- Read lines 42-60 of `plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs`
  to confirm field names.
- Run verification step 4 (BACKLOG.md test) and record the output.
- Search `plugins/*/hooks/*.js` for `statSync` to confirm or deny prior usage.

Revision rule:

If `tool_input.path` is not the correct field name, correct the field access in the implementation
and state the actual field name in the handoff. If `fs.statSync` is found to be inappropriate
(e.g., hook runs before the path exists on disk), fall back to the heuristic-only approach (Option
B from section 5.2) and document why.

---

### T3: Fix SKILL.md Frontmatter and Linting

```yaml
---
task: T3
title: Fix SKILL.md frontmatter user-invocable field and run linting
status: not-started
agent: general-purpose
dependencies: []
priority: 2
complexity: low
accuracy-risk: medium
parallelize-with: [T1, T2]
reason: T3 touches only SKILL.md frontmatter. T1 touches plugin.json and rules/CLAUDE.md. T2 touches the JS hook. No overlap.
handoff: "Report: final frontmatter of SKILL.md (full YAML block), output of `uv run prek run --files`, and rationale for the user-invocable value chosen."
---
```

#### Context

`plugins/orchestrator-discipline/skills/orchestrator-discipline/SKILL.md` (83 lines) has minimal
frontmatter:

```yaml
---
description: "..."
---
```

Two questions are unresolved per the feature context doc:

- Q2: Should `user-invocable: true` be added? (The `commands` field in `plugin.json` points to this
  skills directory, which may register it as a slash command. The `plugin-creator` CLAUDE.md notes
  that a `name:` field causes skills NOT to appear as slash commands — the current SKILL.md has no
  `name:` field, which is correct.)
- Q4: The `commands` field in `plugin.json` declares the same path as `skills`. After T1 resolves
  the `commands` redundancy question, T3 must ensure SKILL.md frontmatter is correct for whichever
  registration mechanism remains.

Per `plugins/plugin-creator/CLAUDE.md` (skill name field bug section, per feature context doc line
282): `user-invocable: true` must be set for the skill to appear in the slash command menu.

#### Objective

Determine the correct value for `user-invocable` in the SKILL.md frontmatter, add the field, and
confirm the file passes `uv run prek run --files`.

#### Required Inputs

- `plugins/orchestrator-discipline/skills/orchestrator-discipline/SKILL.md` — full file to read
  before editing
- `plugins/plugin-creator/CLAUDE.md` — skill name field bug section (determines whether
  `user-invocable: true` is needed and whether `name:` field must be absent)
- `plugins/plugin-creator/skills/claude-plugins-reference-2026/SKILL.md` — SKILL.md frontmatter
  field reference (what fields are valid)
- One working example SKILL.md from another plugin that uses `user-invocable: true` — confirm the
  field name and value format (true vs "true" vs boolean)
- T1 handoff (if available): whether `commands` field was removed from `plugin.json`. If `commands`
  was removed, `user-invocable: true` is the only registration mechanism — it becomes required, not
  optional.

#### Requirements

1. Read the full current SKILL.md before editing.
2. Read `plugins/plugin-creator/CLAUDE.md` to find the definitive answer on whether
   `user-invocable: true` is required for slash command registration when `name:` is absent.
3. Find one existing SKILL.md in this repo that uses `user-invocable: true` and confirm the exact
   YAML syntax (field name, value format).
4. Add `user-invocable: true` to the frontmatter if confirmed required. Place it after
   `description:` and before any other fields, consistent with the example found in requirement 3.
5. Do NOT add a `name:` field — the feature context doc (line 282) explicitly states the absence
   of `name:` is correct.
6. Run `uv run prek run --files plugins/orchestrator-discipline/skills/orchestrator-discipline/SKILL.md`
   and resolve any linting findings.

#### Constraints

- Do NOT change the body content of SKILL.md — only frontmatter.
- Do NOT add `name:` to the frontmatter.
- Do NOT change the `description:` value.
- If T1 has not completed, proceed with `user-invocable: true` as the default choice (the feature
  context doc states manual invocation is a desired use case per scenario 6, line 136-139).
- Do NOT run `git commit`.

#### Expected Outputs

- `plugins/orchestrator-discipline/skills/orchestrator-discipline/SKILL.md` — modified frontmatter
  with `user-invocable: true` added (and any linting fixes applied)

#### Acceptance Criteria

1. SKILL.md frontmatter contains `user-invocable: true` (exact YAML boolean, not string `"true"`).
2. SKILL.md frontmatter does NOT contain a `name:` field.
3. `uv run prek run --files plugins/orchestrator-discipline/skills/orchestrator-discipline/SKILL.md`
   exits 0 with no errors.
4. SKILL.md body content (lines after the closing `---`) is unchanged from the original.

#### Verification Steps

1. Read the modified SKILL.md and print the frontmatter block (lines 1 to closing `---`). Confirm
   `user-invocable: true` is present and `name:` is absent.
2. `uv run prek run --files plugins/orchestrator-discipline/skills/orchestrator-discipline/SKILL.md`
   — record full output; confirm exit 0.
3. `grep -n "user-invocable" plugins/orchestrator-discipline/skills/orchestrator-discipline/SKILL.md`
   — confirm line number and value.
4. `grep -n "^name:" plugins/orchestrator-discipline/skills/orchestrator-discipline/SKILL.md` —
   confirm no output (field absent).

#### CoVe Checks

Key claims to verify:

- Claim: `user-invocable: true` (YAML boolean) is the correct field name and value format for
  slash command registration in plugin SKILL.md files.
- Claim: The absence of a `name:` field is intentional and correct (not a missing required field).

Verification questions:

1. Does any existing SKILL.md in `plugins/*/skills/` use `user-invocable: true`? Find one and read
   its frontmatter to confirm exact syntax. (`Grep(pattern="user-invocable", path="plugins/")`)
2. Does `plugins/plugin-creator/CLAUDE.md` or the plugin-creator skill reference explicitly state
   that `name:` in SKILL.md frontmatter prevents slash command registration?
3. Does `plugins/plugin-creator/skills/claude-plugins-reference-2026/SKILL.md` list `user-invocable`
   as a valid frontmatter field for skills?

Evidence to collect:

- Run `Grep(pattern="user-invocable", path="plugins/")` and read one matching file to confirm
  field syntax.
- Read the relevant section of `plugins/plugin-creator/CLAUDE.md` and cite the exact line.

Revision rule:

If `user-invocable` is not found in any existing SKILL.md, the field name may be wrong. Check the
plugin-creator reference for the correct field name before writing. State the finding and actual
field name in the handoff.

---

## SYNC CHECKPOINT 1: T1 + T2 + T3 Convergence

Convergence point: All three of T1, T2, T3 must complete before T4 starts.

Quality gates:

- T1: `claude plugin validate plugins/orchestrator-discipline/` exits 0
- T1: `plugin.json` contains no `rules` key; `rules/CLAUDE.md` content preserved in verified
  destination
- T2: All 5 verification step commands for T2 produce expected outputs
- T2: `node --check` passes on the modified hook file
- T3: `user-invocable: true` present in SKILL.md frontmatter, `name:` absent
- T3: `uv run prek run --files` passes on SKILL.md

Reflection questions:

- Did T1 remove the `commands` field? If so, does T3's frontmatter still correctly register the
  skill?
- Did T2 introduce any new import (`require`) statements? If so, does T1's updated `plugin.json`
  need any adjustment?
- Are there any emergent conflicts between the three parallel fixes?

Proceed to T4 only after all T1 + T2 + T3 acceptance criteria are confirmed met.

---

### T4: Full Validation Suite

```yaml
---
task: T4
title: Full validation suite — end-to-end verification of all orchestrator-discipline fixes
status: not-started
agent: general-purpose
dependencies: [T1, T2, T3]
priority: 1
complexity: medium
accuracy-risk: low
parallelize-with: []
reason: T4 is the convergence task; it depends on all three prior tasks completing. No parallelism possible.
handoff: "Report: output of each validation command (copy-pasted, not summarized), pass/fail status for each acceptance criterion, any pre-existing issues found in plugin_validator.py output, and the documented test procedure for hook behavior."
---
```

#### Context

T1, T2, and T3 have each addressed one gap in the orchestrator-discipline plugin:

- T1: Removed invalid `rules` field from `plugin.json`
- T2: Extended Grep directory detection in `pre-tool-orchestrator-read-warning.cjs`
- T3: Added `user-invocable: true` to SKILL.md frontmatter

T4 runs the complete validation suite across all plugin files to confirm the fixes integrate
correctly and no regressions exist.

Two validators exist for this repo:

- `claude plugin validate plugins/orchestrator-discipline/` — official Claude Code CLI validator
- `./plugins/plugin-creator/scripts/plugin_validator.py plugins/orchestrator-discipline/` — repo
  internal validator script

Both must be run. Failures in `plugin_validator.py` must be documented; pre-existing failures
(present before T1/T2/T3 changes) should be noted as pre-existing and added to the backlog if not
fixed.

#### Objective

Confirm that all fixes from T1, T2, and T3 integrate correctly, both validators pass (or known
pre-existing failures are documented), all modified files pass `prek` linting, and the hook
behavior is verified end-to-end via the five test commands defined in T2's verification steps.

#### Required Inputs

- All modified plugin files (outputs from T1, T2, T3)
- T1, T2, T3 handoff reports (review before starting)
- `plugins/plugin-creator/scripts/plugin_validator.py` — internal validator
- `uv run prek run` — linting tool

#### Requirements

1. Verify T1, T2, and T3 handoff reports confirm all their acceptance criteria are met before
   proceeding.
2. Run `claude plugin validate plugins/orchestrator-discipline/` and capture full output.
3. Run `./plugins/plugin-creator/scripts/plugin_validator.py plugins/orchestrator-discipline/` and
   capture full output.
4. Collect the list of all files modified across T1, T2, T3 and run
   `uv run prek run --files <each-modified-file>` on each one.
5. Re-run the five hook behavior test commands from T2's verification steps to confirm hook
   behavior end-to-end (not just the modified logic in isolation).
6. Document the test procedure for hook behavior as a short code block that can be re-run by a
   human reviewer.
7. If `plugin_validator.py` reports errors not related to the T1/T2/T3 changes (pre-existing
   issues), document each one with file and line reference, and add a backlog item for each.

#### Constraints

- Do NOT fix issues discovered in step 7 in this task — document them only. Scope boundaries apply.
- Do NOT run `git commit` or `git add`.
- Do NOT modify any files unless a T1/T2/T3 fix was not correctly applied (in which case, apply
  only the missing part of the fix with a note in the handoff).

#### Expected Outputs

- A validation report (written to stdout / handoff) containing:
  - Full output of `claude plugin validate` (copy-pasted verbatim)
  - Full output of `plugin_validator.py` (copy-pasted verbatim)
  - Full output of `prek` linting for each modified file
  - Pass/fail result for each acceptance criterion
  - Hook behavior test results (all five commands and their outputs)
  - List of any pre-existing issues found (if any)

#### Acceptance Criteria

1. `claude plugin validate plugins/orchestrator-discipline/` exits 0 with 0 errors.
2. `plugin_validator.py plugins/orchestrator-discipline/` exits 0, OR any non-zero exit is
   explained by pre-existing issues documented with file:line citations.
3. `uv run prek run --files` passes (exit 0) for every file modified in T1, T2, and T3.
4. Hook test — directory path: `echo '{"tool_name":"Grep","tool_input":{"path":"src/","pattern":"foo"}}' | node plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs` produces JSON containing `"additionalContext"`.
5. Hook test — source file path: `echo '{"tool_name":"Grep","tool_input":{"path":"plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs","pattern":"foo"}}' | node plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs` produces JSON containing `"additionalContext"`.
6. Hook test — markdown path: `echo '{"tool_name":"Grep","tool_input":{"path":"BACKLOG.md","pattern":"foo"}}' | node plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs` produces `{}`.
7. Hook test — empty path: `echo '{"tool_name":"Grep","tool_input":{"path":"","pattern":"foo"}}' | node plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs` produces `{}` and exits without error.
8. SKILL.md frontmatter contains `user-invocable: true` and no `name:` field (confirmed by grep).
9. `plugin.json` contains no `rules` key (confirmed by grep or cat).
10. `rules/CLAUDE.md` content (the "no exemption categories" language) is present in the
    destination file chosen in T1.
11. Hook test — extensionless directory path: `echo '{"tool_name":"Grep","tool_input":{"path":"plugins/orchestrator-discipline/hooks","pattern":"def "}}' | node plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs` produces JSON containing `"additionalContext"` (tests the `looksLikeDirectory` heuristic for paths with no file extension and no trailing slash).

#### Verification Steps

1. `claude plugin validate plugins/orchestrator-discipline/` — copy full output to handoff.
2. `./plugins/plugin-creator/scripts/plugin_validator.py plugins/orchestrator-discipline/` — copy full output to handoff.
3. For each modified file path from T1/T2/T3 handoffs: `uv run prek run --files <path>` — copy output for each.
4. Run all five echo-pipe test commands from acceptance criteria 4–7 above and copy outputs.
5. `grep -n "rules" plugins/orchestrator-discipline/.claude-plugin/plugin.json` — confirm no output.
6. `grep -n "user-invocable" plugins/orchestrator-discipline/skills/orchestrator-discipline/SKILL.md` — confirm `user-invocable: true` on output.
7. `grep -n "no exemption categories" <destination-file-from-T1>` — confirm content present.
8. `echo '{"tool_name":"Grep","tool_input":{"path":"plugins/orchestrator-discipline/hooks","pattern":"def "}}' | node plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs` — output must contain `additionalContext` (extensionless directory path coverage).

---

## SYNC CHECKPOINT 2: Final Review

Convergence point: T4 must complete and all acceptance criteria pass.

Quality gates:

- `claude plugin validate` exits 0
- `plugin_validator.py` exits 0 or pre-existing issues are fully documented
- All prek linting passes
- All five hook behavior tests produce expected outputs
- `rules/CLAUDE.md` content preserved in verified destination
- No unaddressed issues introduced by T1/T2/T3

If any quality gate fails, return the failing task to `not-started` and re-execute with the
specific failure mode documented.

Ready for human review and commit when all gates pass.
