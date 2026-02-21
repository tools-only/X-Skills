# Architecture Spec: Validate and Fix orchestrator-discipline Plugin

**Date**: 2026-02-20
**Type**: Validation and fix (not new feature)
**Scope**: Four targeted changes to make `claude plugin validate` pass, close the Grep directory gap, confirm rules delivery, and resolve SKILL.md registration.

---

## Decision Record

### Decision 1: plugin.json `rules` field

**Chosen option**: Move `rules/CLAUDE.md` to `plugins/orchestrator-discipline/CLAUDE.md` (plugin root), then remove the `"rules"` field from plugin.json and delete the `rules/` directory.

**Rationale**:

The official plugin.json schema (per `plugins/plugin-creator/skills/claude-plugins-reference-2026/SKILL.md:32-52`) lists valid component path fields as: `commands`, `agents`, `skills`, `hooks`, `mcpServers`, `outputStyles`, `lspServers`. The field `rules` is absent. The validator error is `Unrecognized key: "rules"` — direct consequence of this absence.

No other plugin in this repo uses a `rules/` directory or a `rules` field in plugin.json. Two other plugins (`plugins/development-harness/CLAUDE.md`, `plugins/plugin-creator/CLAUDE.md`) use a root-level `CLAUDE.md` at the plugin directory. This is a confirmed, working pattern for delivering behavioral context into sessions where the plugin is loaded.

Option A (remove `"rules"` and leave `rules/CLAUDE.md` undeclared) does not deliver the content — there is no evidence Claude Code auto-discovers `rules/` subdirectories in plugins.

Option B (merge into SKILL.md) conflates two different content types: SKILL.md is activated guidance loaded on demand, while the rules are session-level behavioral constraints that should be active passively whenever the plugin is loaded. Merging them means the constraints only apply when a user explicitly invokes the skill.

Option C (root-level `CLAUDE.md`) matches the confirmed pattern, delivers content passively at plugin load time, and removes the invalid schema field.

**Interface changes**:

`plugins/orchestrator-discipline/.claude-plugin/plugin.json`:

- Remove the `"rules": ["./rules"]` line entirely
- Remove the `"commands": ["./skills/orchestrator-discipline"]` line (see Decision 4 below)

`plugins/orchestrator-discipline/CLAUDE.md`:

- Create this file by moving the full content of `plugins/orchestrator-discipline/rules/CLAUDE.md` into it
- Do not alter the content — the behavioral constraints are correct as written

`plugins/orchestrator-discipline/rules/CLAUDE.md`:

- Delete this file after content is confirmed moved to the plugin root `CLAUDE.md`

`plugins/orchestrator-discipline/rules/` directory:

- Delete the now-empty directory

---

### Decision 2: Grep directory path coverage in hook

**Chosen option**: Add a directory detection branch to the Grep path evaluation in `pre-tool-orchestrator-read-warning.cjs`. Use `fs.statSync(targetPath).isDirectory()` as the primary check, with `looksLikeDirectory()` heuristic as fallback. Fire the warning for all Grep calls where the path resolves to a directory — no source-name allowlist.

**Rationale**:

The options from the feature context were:

- A: Fire on all Grep directory paths
- B: Fire on directory AND glob/pattern suggests source files
- C: Fire except for an allowed list of non-source directories
- D: Apply extension check to the Grep pattern parameter

Option B and D require reasoning about the Grep pattern (e.g., `*.py` implies source code, `*.md` implies docs). This adds complexity and can still miss extension-agnostic patterns like `class.*Service` or `TODO`. Pattern-content heuristics are fragile.

Option C (allowlist) creates maintenance burden: `plan/`, `docs/`, `references/` today; next month a new directory. Allowlists drift.

Option A fires on all directory Grep calls. The concern raised in the feature context is noise from legitimate orchestrator operations like `Grep(path="plan/")`. This concern is partially mitigated by: (a) the hook is non-blocking — it injects `additionalContext`, not a block, and (b) the warning text includes the falsifiable test ("Will I Edit or Write this file?") which is the correct self-check even for plan files. An orchestrator Grep on `plan/` to find a task ID is a legitimate read-only operation, but the hook firing is a nudge, not a prohibition. The noise cost is acceptable given the benefit of catching directory-scoped source investigations.

`fs.statSync` is the reliable approach when the path exists. Grep targets always exist (if they do not, Grep fails before the hook fires). The `looksLikeDirectory()` heuristic (ends with `/`, or last segment has no `.`) is the fallback for edge cases such as relative paths that resolve outside the hook's working directory.

**Exact logic change**:

Location: `plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs`

Add `require('node:fs')` at the top of the file alongside existing requires (the file currently uses no `fs` import).

Add two new functions after the existing `isSourceOrConfigFile` and `isTestPath` functions:

```text
isDirectory(targetPath):
  try { return fs.statSync(targetPath).isDirectory(); } catch { return false; }

looksLikeDirectory(targetPath):
  if (!targetPath) return false
  if targetPath ends with '/' or '\\': return true
  lastSegment = targetPath split on path separators, last element
  return lastSegment does not contain '.'
```

Change the Grep branch of the existing path-evaluation logic from:

```text
if toolName is Grep:
  targetPath = toolInput.path || ''
  shouldWarn = isSourceOrConfigFile(targetPath) || isTestPath(targetPath)
```

to:

```text
if toolName is Grep:
  targetPath = toolInput.path || ''
  shouldWarn = isSourceOrConfigFile(targetPath)
              || isTestPath(targetPath)
              || isDirectory(targetPath)
              || looksLikeDirectory(targetPath)
```

The Read branch is unchanged.

**additionalContext output for directory Grep**:

The hook already produces a single `<orchestrator-read-warning>` block regardless of how the warning was triggered. No changes to the output format are needed. The existing warning text — which references "Grep" explicitly alongside "Read" — is already accurate for the directory case.

**Risk**: The `fs` module is already available in Node.js CommonJS context (`require('node:fs')`). No new dependencies. `statSync` throws if the path does not exist; the try/catch handles this. False negatives (missed directories) are bounded by `looksLikeDirectory` fallback.

---

### Decision 3: SKILL.md `user-invocable` field

**Chosen option**: Add `user-invocable: true` to the SKILL.md frontmatter.

**Rationale**:

The feature context (line 283) cites the skill name field bug: "Plugin skills with an explicit `name:` field in SKILL.md frontmatter do NOT appear as slash commands. The SKILL.md at `skills/orchestrator-discipline/SKILL.md` currently has no `name:` field — this is correct. `user-invocable: true` must be set for the skill to appear in the slash command menu."

The plugin currently sets both `"skills"` and `"commands"` in plugin.json pointing to the same path. Decision 4 removes `"commands"`. After that removal, `user-invocable: true` in SKILL.md is the mechanism that registers the skill as a slash command.

Whether passive activation (loading skill context without explicit invocation) is the only use case: the plugin's hooks run automatically. But the skill contains the discipline guidance text — the rationale, the falsifiable test, the anti-patterns. A user who wants to review the full discipline framework should be able to invoke `/orchestrator-discipline:orchestrator-discipline` directly. This is a valid deliberate activation scenario (use scenario 6 in the feature context).

**Interface change**:

`plugins/orchestrator-discipline/skills/orchestrator-discipline/SKILL.md` frontmatter:

Current:

```yaml
---
description: "..."
---
```

After:

```yaml
---
description: "..."
user-invocable: true
---
```

No `name:` field is added (adding `name:` would suppress slash command registration per the bug note).

---

### Decision 4: Remove `commands` field from plugin.json

**Chosen option**: Remove `"commands": ["./skills/orchestrator-discipline"]` from plugin.json.

**Rationale**:

The official plugin.json schema and the `agent-orchestration` plugin example both use `skills` as the current field. The `commands` field is documented as legacy. Double-registration of the same path under both `skills` and `commands` is undocumented behavior. The codebase patterns analysis (section 3.2) shows the current plugin.json has both fields pointing to identical paths.

After removing `commands`, the skill registration path is: `"skills": ["./skills/orchestrator-discipline"]` in plugin.json, plus `user-invocable: true` in SKILL.md frontmatter (Decision 3). This is the canonical, non-legacy registration pattern.

**Interface change**:

`plugins/orchestrator-discipline/.claude-plugin/plugin.json` after both Decision 1 and Decision 4:

```json
{
  "name": "orchestrator-discipline",
  "skills": ["./skills/orchestrator-discipline"],
  "hooks": "./hooks.json"
}
```

---

## Exact Interface Changes Summary

### File: `plugins/orchestrator-discipline/.claude-plugin/plugin.json`

Remove fields: `"rules"`, `"commands"`

Resulting structure:

```json
{
  "name": "orchestrator-discipline",
  "skills": ["./skills/orchestrator-discipline"],
  "hooks": "./hooks.json"
}
```

### File: `plugins/orchestrator-discipline/CLAUDE.md` (new)

Content: full content of `plugins/orchestrator-discipline/rules/CLAUDE.md`, verbatim, no modifications.

### File: `plugins/orchestrator-discipline/rules/CLAUDE.md` (delete)

After content is confirmed moved to plugin root `CLAUDE.md`.

### File: `plugins/orchestrator-discipline/hooks/pre-tool-orchestrator-read-warning.cjs`

Three changes:

1. Add `const fs = require('node:fs');` to the require block at the top of the file.

2. Add two new helper functions after the existing `isTestPath` function:

   - `isDirectory(targetPath)` — uses `fs.statSync(targetPath).isDirectory()` in a try/catch returning false on error
   - `looksLikeDirectory(targetPath)` — returns true if path ends with `/` or `\`, or if the last path segment contains no `.`

3. In the Grep branch where `shouldWarn` is evaluated: extend the condition to also return true when `isDirectory(targetPath)` or `looksLikeDirectory(targetPath)` is true.

No changes to: output format, Read branch logic, the `additionalContext` XML content, exit codes, or error handling.

### File: `plugins/orchestrator-discipline/skills/orchestrator-discipline/SKILL.md`

Add `user-invocable: true` to YAML frontmatter. No `name:` field.

---

## Acceptance Criteria

| # | Criterion | Verification Method |
|---|-----------|---------------------|
| 1 | `claude plugin validate plugins/orchestrator-discipline/` exits with no errors | Run command, observe zero error lines |
| 2 | Hook fires when `Read(file_path="src/cli/commands.py")` — source file path | Simulate hook stdin, confirm `additionalContext` in output |
| 3 | Hook fires when `Grep(pattern="class", path="src/")` — directory path | Simulate hook stdin with directory path, confirm `additionalContext` in output |
| 4 | Hook fires when `Grep(pattern="def ", path="plugins/orchestrator-discipline/hooks")` — extensionless directory path without trailing slash | Confirm `looksLikeDirectory` or `isDirectory` triggers |
| 5 | Hook does NOT fire for `Grep(pattern="ADR", path="plan/feature-context.md")` — single `.md` file | Confirm empty `{}` output |
| 6 | Hook does NOT fire for `Read(file_path="BACKLOG.md")` — doc file not matching source extensions | Confirm empty `{}` output |
| 7 | `rules/CLAUDE.md` content exists verbatim in `plugins/orchestrator-discipline/CLAUDE.md` | Diff old vs new file, no content loss |
| 8 | `plugins/orchestrator-discipline/rules/` directory no longer exists | `ls plugins/orchestrator-discipline/` shows no `rules/` entry |
| 9 | `plugin.json` contains no `"rules"` key and no `"commands"` key | `grep -c rules .claude-plugin/plugin.json` returns 0 |
| 10 | SKILL.md frontmatter includes `user-invocable: true` and no `name:` field | Read frontmatter, confirm both conditions |
| 11 | `uv run prek run --files <changed files>` passes with no errors | Run prek on all modified files |
| 12 | Commit message follows `fix(orchestrator-discipline): <description>` format | `git log --oneline -1` |

---

## Risk Notes

### Risk 1: Rules content delivery timing

Moving `rules/CLAUDE.md` to a plugin root `CLAUDE.md` relies on Claude Code loading that file when the plugin is activated. The pattern is observed in `plugins/development-harness/CLAUDE.md` and `plugins/plugin-creator/CLAUDE.md`, but the exact load timing (at install, at session start, or on demand) is not documented in the codebase analysis. The hooks deliver equivalent enforcement on every tool call regardless of whether the CLAUDE.md loads. If CLAUDE.md load timing is uncertain, the behavioral constraints are still enforced structurally via hooks. The CLAUDE.md provides the rationale text; hooks provide the enforcement.

**Mitigation**: This is a net improvement over the current state where the `rules` field is unrecognized and may never have loaded at all.

### Risk 2: Grep directory hook noise

The `isDirectory` + `looksLikeDirectory` approach fires the warning for any directory path in a Grep call. An orchestrator Grep on `plan/` or `docs/` will receive the read-warning. The warning is non-blocking and its text is the correct self-check even in those cases. Over time, if the noise is observed to degrade signal, an allowlist can be added as a follow-on patch. The spec does not include an allowlist to keep the initial fix minimal.

### Risk 3: `commands` removal and double-registration

Removing `"commands"` from plugin.json when `"skills"` already registers the same path is safe per the schema docs. The only risk is if Claude Code uses `commands` as the mechanism for slash command registration and ignores `skills` for that purpose. Decision 3 (`user-invocable: true` in SKILL.md) is the safeguard: it explicitly opts the skill into slash command registration via the skills mechanism, independent of the `commands` field.

### Risk 4: Auto-sync version bump collision

The `auto_sync_manifests.py` pre-commit hook bumps plugin version on file changes. The implementer must not manually edit version fields in plugin.json or marketplace.json. Allow the pre-commit hook to handle versioning. If multiple files in the plugin change in a single commit, auto-sync runs once and produces one version bump.

### Risk 5: `fs.statSync` path resolution

The hook executes in the context of the Claude Code session's working directory. `fs.statSync(targetPath)` resolves relative paths against the hook process's working directory. If the Grep target path is relative (e.g., `src/`), the stat call resolves it relative to the session's cwd. This should be the repo root in normal Claude Code operation. If cwd is unexpected, `statSync` throws and the catch block returns false — `looksLikeDirectory` heuristic applies as fallback. No data is lost; at worst the warning does not fire for that specific path format.
