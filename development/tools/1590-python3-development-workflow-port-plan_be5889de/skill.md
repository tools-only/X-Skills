# Python3-Development Plugin: SAM Workflow Port Plan (from gitlab-runner-management)

This document is a **durable execution plan** for porting the “add-new-feature → implement-feature → start-task → complete-implementation” workflow (and its supporting agents/scripts) from `~/repos/gitlab-runner-management` into this repository’s `plugins/python3-development/` plugin, **as skills** (commands deprecated).

The goal is that you can **stop at any checkpoint** and resume later without re-deriving context.

---

## Goals

- Port the **artifact-driven, gated workflow** from `gitlab-runner-management` into `plugins/python3-development` using **skills + agents + hooks**.
- Make `plugins/python3-development` **self-contained** (no references to `~/.claude/...` files that aren’t in this repo).
- Remove/avoid reliance on `plugins/python3-development/commands/` (or keep as legacy reference only, but not required for execution).

---

## Inputs (Source of Truth)

### Reference repo (ground truth)

- `~/repos/gitlab-runner-management/.claude/commands/`
  - `add-new-feature.md`
  - `implement-feature.md`
  - `start-task.md`
  - `complete-implementation.md`
- `~/repos/gitlab-runner-management/.claude/agents/`
  - `feature-researcher.md`
  - `codebase-analyzer.md`
  - `plan-validator.md`
  - `context-gathering.md`
  - `feature-verifier.md`
  - (plus any “complete-implementation phase” agents it calls: `integration-checker.md`, `doc-drift-auditor.md`, `service-documentation.md`, `context-refinement.md`, `code-reviewer.md`)
- `~/repos/gitlab-runner-management/.claude/skills/implementation-manager/`
  - `SKILL.md`
  - `scripts/implementation_manager.py`
  - `scripts/task_status_hook.py`

### Target plugin (where changes land)

- `./plugins/python3-development/.claude-plugin/plugin.json`
- `./plugins/python3-development/skills/`
- `./plugins/python3-development/agents/`

---

## Checkpoint 0 — Safety + State Capture (5 minutes)

- [ ] Run `git status` in this repo and record:
  - current branch
  - untracked/modified files
- [ ] Create a new working branch (recommended) for the port
- [ ] Decide whether to commit unrelated untracked files before starting

**Resume point:** If interrupted before Checkpoint 1, just rerun `git status` and continue.

---

## Checkpoint 1 — Inventory: map “workflow components” to plugin equivalents

### 1A) Define the new skill set (user-invocable)

Create these skill directories under `plugins/python3-development/skills/development/`:

- `add-new-feature/` (user-invocable)
- `implement-feature/` (user-invocable)
- `start-task/` (user-invocable)
- `complete-implementation/` (user-invocable)

Each must have `SKILL.md` with:

- `name`: kebab-case (matches directory)
- `description`: include trigger keywords (“add feature”, “plan feature”, “implement feature”, etc.)
- `argument-hint`: modeled after the reference commands
- `user-invocable: true`
- Hooks (where needed), per Checkpoint 4

### 1B) Define support skills (non-user-invocable)

- `skills/implementation-manager/` (user-invocable should be `false`)
  - Provides `implementation_manager.py` and `task_status_hook.py`
  - Provides documentation for JSON outputs + task file schema

### 1C) Define required agents (plugin-bundled)

Port agents from `gitlab-runner-management/.claude/agents/` into `plugins/python3-development/agents/`:

- Required for planning pipeline:
  - `feature-researcher.md`
  - `codebase-analyzer.md`
  - `plan-validator.md`
  - `context-gathering.md`
- Required for completion/verification pipeline:
  - `feature-verifier.md`
  - `integration-checker.md`
  - `doc-drift-auditor.md`
  - `service-documentation.md`
  - `context-refinement.md`
  - `code-reviewer.md` (or map to existing `python-code-reviewer.md` if the prompts are compatible)

**Resume point:** If interrupted after Checkpoint 1, you should have a complete “bill of materials” list.

---

## Checkpoint 2 — Vendor files into this repo (copy, don’t rewrite yet)

### 2A) Vendor the implementation-manager skill + scripts

Copy into the plugin:

- From: `~/repos/gitlab-runner-management/.claude/skills/implementation-manager/SKILL.md`
  - To: `plugins/python3-development/skills/implementation-manager/SKILL.md`
- From: `~/repos/gitlab-runner-management/.claude/skills/implementation-manager/scripts/implementation_manager.py`
  - To: `plugins/python3-development/skills/implementation-manager/scripts/implementation_manager.py`
- From: `~/repos/gitlab-runner-management/.claude/skills/implementation-manager/scripts/task_status_hook.py`
  - To: `plugins/python3-development/skills/implementation-manager/scripts/task_status_hook.py`

### 2B) Vendor the phase agents

Copy the required `.md` agent files into `plugins/python3-development/agents/`.

### 2C) Vendor any templates referenced by the python3 plugin that are currently external/missing

These are known gaps in `plugins/python3-development` today:

- `~/.claude/templates/feature-task-template.md` (referenced by old `create-feature-task` command docs)
- `templates/test-checklist.md` (referenced via `@include` in old command docs)
- `~/.claude/agents/python-cli-demo.py` (referenced by multiple python3 docs/agents)
- `~/.claude/tools/validate_pep723.py` (referenced by python3 tool registry docs)

Pick **one canonical location** inside `plugins/python3-development/` (recommendation):

- Templates:
  - `plugins/python3-development/skills/development/add-new-feature/templates/`
  - `plugins/python3-development/skills/testing/comprehensive-test-review/templates/`
- Example code/assets:
  - `plugins/python3-development/skills/python3-development/assets/`
- Helper scripts:
  - `plugins/python3-development/scripts/` (or `skills/python3-development/scripts/`)

**Resume point:** After Checkpoint 2, everything exists locally even if not yet wired.

---

## Checkpoint 3 — Convert reference “commands” into plugin “skills”

For each of these reference command files:

- `gitlab-runner-management/.claude/commands/add-new-feature.md`
- `gitlab-runner-management/.claude/commands/implement-feature.md`
- `gitlab-runner-management/.claude/commands/start-task.md`
- `gitlab-runner-management/.claude/commands/complete-implementation.md`

Create a corresponding `SKILL.md` in the target directories from Checkpoint 1A.

**Critical adaptations when converting:**

- Replace references like:

  - `$CLAUDE_PROJECT_DIR/.claude/skills/implementation-manager/...`
  - `.claude/agents/...`

  with plugin-bundled equivalents (hooks/scripts should reference `${CLAUDE_PLUGIN_ROOT}` where appropriate).

- Keep `$ARGUMENTS` usage intact.
- Keep structured phase outputs and DONE/BLOCKED signaling intact.

**Resume point:** After Checkpoint 3, user can invoke `/add-new-feature`-equivalent skills from the plugin.

---

## Checkpoint 4 — Hook wiring (task status + activity timestamps)

Implement hook wiring equivalent to the reference repo:

- **On `implement-feature` skill**:
  - `SubagentStop` hook runs `task_status_hook.py`
- **On `start-task` skill**:
  - `PostToolUse` hook runs the hook script that updates LastActivity

Make sure the hook command points to the vendored script path inside the plugin. Prefer `${CLAUDE_PLUGIN_ROOT}` for plugin-relative execution.

**Resume point:** After Checkpoint 4, the status automation works again (core “don’t get stuck” feature).

---

## Checkpoint 5 — Update plugin manifest (`plugin.json`)

Update `plugins/python3-development/.claude-plugin/plugin.json`:

- Add the new skills:
  - `./skills/development/add-new-feature`
  - `./skills/development/implement-feature`
  - `./skills/development/start-task`
  - `./skills/development/complete-implementation`
  - `./skills/implementation-manager`
- Add the new agents (individual `./agents/*.md` paths)
- Decide policy for `commands/`:
  - Preferred: do not register `commands` at all (skills-only plugin)

**Resume point:** After Checkpoint 5, plugin installation will expose the new workflow.

---

## Checkpoint 6 — Remove/neutralize legacy `commands/` dependencies in this plugin

The existing `plugins/python3-development/commands/` currently contains `@include templates/...` references that don’t exist in this plugin.

Choose one:

- **Option A (preferred)**: remove `commands/` from “required execution path” entirely
  - Keep the folder only as human documentation (no includes; no external references)
- **Option B**: fully vendor the referenced templates and keep commands working, but treat them as legacy

**Resume point:** After Checkpoint 6, the plugin has no dangling references.

---

## Checkpoint 7 — Verification (don’t skip)

- [ ] Run repo linting on all touched files:
  - `uv run prek run --files <changed-files...>`
- [ ] Validate `plugin.json` is valid JSON
- [ ] Validate the skill frontmatter loads correctly (name/description/argument-hint formatting)
- [ ] Spot check that no `~/.claude/...` file references remain unless intentionally kept as “optional external reference”

**Resume point:** After Checkpoint 7, you should be safe to ship.

---

## Definition of Done (for this port)

- [ ] The planning/execution/verification workflow exists as **skills** in `plugins/python3-development`
- [ ] All referenced scripts/templates/examples are **vendored into this repo**
- [ ] No required step depends on `~/.claude/...`
- [ ] Hooks update task status/timestamps the same way as the reference workflow
- [ ] `uv run prek run --files ...` passes for modified files
