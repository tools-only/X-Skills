# Architecture Pattern Study

## get-shit-done Pattern

### File Tree

```text
get-shit-done/
├── agents/
│   ├── gsd-codebase-mapper.md
│   ├── gsd-debugger.md
│   ├── gsd-executor.md
│   ├── gsd-integration-checker.md
│   ├── gsd-phase-researcher.md
│   ├── gsd-plan-checker.md
│   ├── gsd-planner.md
│   ├── gsd-project-researcher.md
│   ├── gsd-research-synthesizer.md
│   ├── gsd-roadmapper.md
│   └── gsd-verifier.md
├── commands/
│   └── gsd/
│       ├── add-phase.md
│       ├── add-todo.md
│       ├── audit-milestone.md
│       ├── check-todos.md
│       ├── cleanup.md
│       ├── complete-milestone.md
│       ├── debug.md
│       ├── discuss-phase.md
│       ├── execute-phase.md
│       ├── health.md
│       ├── help.md
│       ├── insert-phase.md
│       ├── list-phase-assumptions.md
│       ├── map-codebase.md
│       ├── new-milestone.md
│       ├── new-project.md
│       ├── pause-work.md
│       ├── plan-milestone-gaps.md
│       ├── plan-phase.md
│       ├── progress.md
│       ├── quick.md
│       ├── remove-phase.md
│       ├── research-phase.md
│       ├── resume-work.md
│       ├── set-profile.md
│       ├── settings.md
│       ├── update.md
│       └── verify-work.md
├── get-shit-done/
│   ├── bin/
│   │   ├── gsd-tools.cjs         (Node.js utility — state, roadmap, commit ops)
│   │   └── gsd-tools.test.cjs
│   ├── references/
│   │   ├── checkpoints.md
│   │   ├── continuation-format.md
│   │   ├── decimal-phase-calculation.md
│   │   ├── git-integration.md
│   │   ├── git-planning-commit.md
│   │   ├── model-profile-resolution.md
│   │   ├── model-profiles.md
│   │   ├── phase-argument-parsing.md
│   │   ├── planning-config.md
│   │   ├── questioning.md
│   │   ├── tdd.md
│   │   ├── ui-brand.md
│   │   └── verification-patterns.md
│   ├── templates/
│   │   ├── codebase/
│   │   │   └── (architecture.md, concerns.md, conventions.md, integrations.md, stack.md, structure.md, testing.md)
│   │   ├── config.json
│   │   ├── context.md
│   │   ├── continue-here.md
│   │   ├── DEBUG.md
│   │   ├── debug-subagent-prompt.md
│   │   ├── discovery.md
│   │   ├── milestone.md
│   │   ├── milestone-archive.md
│   │   ├── phase-prompt.md
│   │   ├── planner-subagent-prompt.md
│   │   ├── project.md
│   │   ├── requirements.md
│   │   ├── research.md
│   │   ├── research-project/
│   │   │   └── (ARCHITECTURE.md, FEATURES.md, PITFALLS.md, STACK.md, SUMMARY.md)
│   │   ├── roadmap.md
│   │   ├── state.md
│   │   ├── summary.md
│   │   ├── summary-complex.md
│   │   ├── summary-minimal.md
│   │   ├── summary-standard.md
│   │   ├── UAT.md
│   │   ├── user-setup.md
│   │   └── verification-report.md
│   └── workflows/
│       ├── add-phase.md
│       ├── audit-milestone.md
│       ├── check-todos.md
│       ├── cleanup.md
│       ├── complete-milestone.md
│       ├── diagnose-issues.md
│       ├── discovery-phase.md
│       ├── discuss-phase.md
│       ├── execute-phase.md
│       ├── execute-plan.md
│       ├── health.md
│       ├── help.md
│       ├── insert-phase.md
│       ├── list-phase-assumptions.md
│       ├── map-codebase.md
│       ├── new-milestone.md
│       ├── new-project.md
│       ├── pause-work.md
│       ├── plan-milestone-gaps.md
│       ├── plan-phase.md
│       ├── progress.md
│       ├── quick.md
│       ├── remove-phase.md
│       ├── research-phase.md
│       ├── resume-project.md
│       ├── set-profile.md
│       ├── settings.md
│       ├── transition.md
│       ├── update.md
│       ├── verify-phase.md
│       └── verify-work.md
├── hooks/
│   ├── gsd-check-update.js
│   └── gsd-statusline.js
├── scripts/
│   └── build-hooks.js
└── package.json
```

### Agent Structure

11 agents, all in `agents/` at repo root. Each is a `.md` file with YAML frontmatter + body. No auto-loading via a `skills` field — they are referenced by commands and spawned via the `Task` tool.

Agent roster with scope:

| Agent | Purpose | Tools |
|---|---|---|
| `gsd-executor` | Executes PLAN.md files atomically, per-task commits, checkpoint handling, SUMMARY.md creation, STATE.md updates | Read, Write, Edit, Bash, Grep, Glob |
| `gsd-planner` | Creates PLAN.md files from phase requirements; gap closure and revision modes | Read, Write, Bash, Glob, Grep, WebFetch, mcp__context7__* |
| `gsd-verifier` | Goal-backward verification of completed phases — checks codebase, not plans | Read, Write, Bash, Grep, Glob |
| `gsd-debugger` | Scientific hypothesis-based bug investigation with persistent debug file state | Read, Write, Edit, Bash, Grep, Glob, WebSearch |
| `gsd-phase-researcher` | Researches technical domain before planning; writes RESEARCH.md consumed by planner | Read, Write, Bash, Grep, Glob, WebSearch, WebFetch, mcp__context7__* |
| `gsd-plan-checker` | Verifies plans WILL achieve goal before execution (goal-backward analysis of plan files, not code) | Read, Bash, Glob, Grep |
| `gsd-roadmapper` | Creates ROADMAP.md and STATE.md from requirements; maps requirements to phases | Read, Write, Bash, Glob, Grep |
| `gsd-codebase-mapper` | Explores codebase and writes structured analysis docs (STACK.md, ARCHITECTURE.md, etc.) for one focus area | Read, Bash, Grep, Glob, Write |
| `gsd-project-researcher` | Researches domain ecosystem before new project; writes STACK.md, FEATURES.md, ARCHITECTURE.md, PITFALLS.md | Read, Write, Bash, Grep, Glob, WebSearch, WebFetch, mcp__context7__* |
| `gsd-research-synthesizer` | Reads 4 parallel researcher outputs; synthesizes into SUMMARY.md; commits all research files | Read, Write, Bash |
| `gsd-integration-checker` | Verifies cross-phase integration and E2E flows; checks connections, not existence | Read, Bash, Grep, Glob |

Each agent has a `color` field and a descriptive `description` that acts as a dispatch hint for Claude (who to use when). Agents declare their spawner explicitly in their `<role>` tag: `Spawned by /gsd:execute-phase orchestrator.`

### Agent + Skill Relationship

No `skills` field in agent frontmatter. There is no auto-loading mechanism. Instead, agents explicitly reference domain knowledge files using `@`-references inside their body:

```text
<execution_context>
@~/.claude/get-shit-done/workflows/execute-phase.md
@~/.claude/get-shit-done/references/ui-brand.md
</execution_context>
```

The `@~/.claude/get-shit-done/references/checkpoints.md` reference inside `gsd-executor` functions as an on-demand reference: `For full automation-first patterns, server lifecycle, CLI handling: See @~/.claude/get-shit-done/references/checkpoints.md`. The agent reads it when it needs it, not automatically.

Templates are referenced similarly: `@~/.claude/get-shit-done/templates/summary.md`

There is NO `skills` field, NO `SKILL.md` pattern, and NO Claude Code marketplace plugin structure in this repo. This is a direct npm-installed system that copies files into `~/.claude/`.

### User-Invocable Skills (Entry Points)

Commands live in `commands/gsd/`. Each command file has YAML frontmatter with `name`, `description`, `argument-hint`, and sometimes `agent` and `allowed-tools`.

Two patterns exist:

**Pattern 1: Thin orchestrator command** (most commands)

The command file is a thin wrapper. It:
1. Declares `allowed-tools` that include `Task`
2. Has an `<execution_context>` block that loads a workflow file via `@~/.claude/get-shit-done/workflows/<name>.md`
3. Has a minimal `<process>` block that says "follow the workflow"
4. Passes `$ARGUMENTS` and live project files (`@.planning/ROADMAP.md`, `@.planning/STATE.md`) as context

Example — `commands/gsd/execute-phase.md`:

```yaml
---
name: gsd:execute-phase
description: Execute all plans in a phase with wave-based parallelization
argument-hint: "<phase-number> [--gaps-only]"
allowed-tools:
  - Read, Write, Edit, Glob, Grep, Bash, Task, TodoWrite, AskUserQuestion
---
<objective>...</objective>
<execution_context>
@~/.claude/get-shit-done/workflows/execute-phase.md
@~/.claude/get-shit-done/references/ui-brand.md
</execution_context>
<context>
Phase: $ARGUMENTS
@.planning/ROADMAP.md
@.planning/STATE.md
</context>
<process>
Execute the execute-phase workflow from @~/.claude/get-shit-done/workflows/execute-phase.md end-to-end.
</process>
```

**Pattern 2: Agent-declaring command** (less common)

The command has an `agent: gsd-planner` field in frontmatter. This tells Claude to hand off to that named agent. Example — `commands/gsd/plan-phase.md`:

```yaml
---
name: gsd:plan-phase
description: Create detailed phase plan (PLAN.md) with verification loop
argument-hint: "[phase] [--auto] [--research] [--skip-research] [--gaps] [--skip-verify]"
agent: gsd-planner
allowed-tools:
  - Read, Write, Bash, Glob, Grep, Task, WebFetch, mcp__context7__*
---
<execution_context>
@~/.claude/get-shit-done/workflows/plan-phase.md
</execution_context>
<process>
Execute the plan-phase workflow from @~/.claude/get-shit-done/workflows/plan-phase.md end-to-end.
</process>
```

### The Handoff Pattern (traced)

Trace for `/gsd:execute-phase 3`:

1. **User types** `/gsd:execute-phase 3`
2. **Claude loads** `commands/gsd/execute-phase.md` — reads the command file into context. The `@`-references cause Claude to also load `workflows/execute-phase.md` and `references/ui-brand.md` and the live project files.
3. **Command body instructs Claude** (the orchestrator) to "Execute the execute-phase workflow end-to-end." The orchestrator IS Claude — it reads the workflow and follows it.
4. **Workflow logic** (in `workflows/execute-phase.md`) tells the orchestrator to: discover plans, analyze dependency waves, then spawn subagents per plan using the Task tool.
5. **For each plan**, orchestrator constructs a prompt using `@~/.claude/get-shit-done/workflows/execute-plan.md` as context and spawns:
   ```
   Task(subagent_type="gsd-executor", prompt=<constructed prompt>)
   ```
6. **`gsd-executor` agent** receives the prompt with the PLAN.md content and `execute-plan.md` workflow loaded. Executes tasks, commits, creates SUMMARY.md, updates STATE.md, returns structured completion message.
7. **Orchestrator** receives the structured return from each agent, handles checkpoints if any, updates routing.

Trace for `/gsd:debug "login fails"`:

1. **User types** `/gsd:debug "login fails"`
2. **Claude loads** `commands/gsd/debug.md` — lightweight command.
3. **Orchestrator** (Claude running the command) gathers symptoms via `AskUserQuestion`, then spawns:
   ```
   Task(subagent_type="gsd-debugger", model=resolved_model, prompt=filled_prompt)
   ```
4. **`gsd-debugger` agent** receives symptoms in `<symptoms>` tags and `symptoms_prefilled: true` mode flag. Runs scientific investigation loop. Returns one of: ROOT CAUSE FOUND, CHECKPOINT REACHED, INVESTIGATION INCONCLUSIVE.
5. **Orchestrator** handles each return type: presents to user, spawns continuation agent if checkpoint, offers options if complete.

Key insight: The WORKFLOW file is the decision tree logic for the orchestrator. The AGENT file is the specialized executor that knows its own domain deeply.

### Skill Categories Found

Three distinct categories, none called "skills" formally:

**1. Commands** (`commands/gsd/*.md`) — User-invocable entry points
- Thin, mostly pass-through
- Declare allowed tools and pass `$ARGUMENTS`
- Load workflow file via `@` reference
- May declare `agent:` for direct handoff

**2. Workflows** (`get-shit-done/workflows/*.md`) — Orchestrator decision trees
- Loaded by commands via `@` references
- Contain the multi-step process logic
- Tell the orchestrator WHAT agents to spawn, WHEN, with WHAT prompt
- Not user-invocable directly — loaded implicitly by commands

**3. References** (`get-shit-done/references/*.md`) — Domain knowledge
- Loaded on-demand by agents via `@` references in their body
- Contain specific technical rules (checkpoint patterns, UI brand, TDD patterns)
- Not agents, not commands — pure reference material

**4. Templates** (`get-shit-done/templates/*.md`) — Output scaffolds
- Agents reference these to produce consistent PLAN.md, SUMMARY.md, VERIFICATION.md output
- Provide the exact structure agents must fill in

**5. Agent files** (`agents/*.md`) — Specialized executors
- Deep domain knowledge embedded directly in the agent body
- Contain the full operational protocol: steps, deviation rules, output formats, success criteria

### Output Contracts

Output formats are defined WITHIN each agent's body. Every agent has a `<structured_returns>` or `<completion_format>` or `<output>` section that defines exactly what it returns to the orchestrator. Example from `gsd-executor`:

```markdown
<completion_format>
## PLAN COMPLETE
**Plan:** {phase}-{plan}
**Tasks:** {completed}/{total}
**SUMMARY:** {path to SUMMARY.md}
**Commits:** {hash}: {message}
**Duration:** {time}
</completion_format>
```

The orchestrator downstream knows what format to expect and pattern-matches on the heading (`## PLAN COMPLETE`, `## CHECKPOINT REACHED`, etc.).

### Routing Mechanism

Routing is encoded in the workflow files, not in a registry. The workflow for `plan-phase.md` explicitly describes when to spawn `gsd-phase-researcher` vs `gsd-plan-checker` vs `gsd-planner` based on flags and state. Routing decision points include:

- Flags from `$ARGUMENTS` (`--research`, `--gaps`, `--skip-verify`)
- State from `gsd-tools.cjs` JSON output (has_research, has_plans, plan_checker_enabled)
- Return values from spawned agents (`## VERIFICATION PASSED` vs `## ISSUES FOUND`)

The orchestrator (Claude) reads the workflow and makes routing decisions inline. There is no separate registry or router file.

---

## claude-code-plugins Pattern

### File Tree

```text
claude-code-plugins/
├── CLAUDE.md
├── CHANGELOG.md
├── README.md
├── .claude-plugin/
│   └── marketplace.json        (registry of all plugins in this repo)
├── rust-developer/
│   ├── .claude-plugin/
│   │   └── plugin.json         (plugin manifest)
│   ├── CLAUDE.md
│   ├── README.md
│   ├── agents/
│   │   └── rust-developer.md
│   ├── commands/
│   │   └── rust-developer.md
│   ├── hooks/
│   │   ├── hooks.json
│   │   └── scripts/
│   │       └── check-rust-patterns.sh
│   └── skills/
│       └── rust-knowledge/
│           ├── SKILL.md
│           ├── assets/
│           │   ├── Cargo.toml.template
│           │   ├── clippy.toml
│           │   └── rustfmt.toml
│           ├── references/
│           │   ├── async-rust.md
│           │   ├── cargo.md
│           │   ├── clippy.md
│           │   ├── error-handling.md
│           │   ├── frameworks.md
│           │   ├── ownership.md
│           │   ├── testing.md
│           │   ├── troubleshooting.md
│           │   └── wasm.md
│           └── scripts/
│               ├── format_lint.sh
│               ├── new_project.sh
│               └── run_tests.sh
├── swift-developer/            (same structure as rust-developer)
├── sql-developer/              (same structure)
└── powershell-developer/       (same structure)
```

Each plugin follows an identical template:
`agents/` + `commands/` + `hooks/` + `skills/{skill-name}/` + `plugin.json`

### Agent Structure

4 plugins, each with exactly 1 agent. All agents have `model: inherit` and `color` in frontmatter.

Agents have rich descriptions with 4 `<example>` blocks in the frontmatter. These examples function as dispatch guidance — they teach Claude WHEN to invoke this agent. Example from `rust-developer`:

```yaml
---
name: rust-developer
description: |
  Use this agent when the user explicitly requests help with Rust development tasks...
  <example>
  Context: User is starting a new Rust project...
  user: "Help me create a new Rust library for handling HTTP requests"
  assistant: "I'll use the rust-developer agent..."
  <commentary>...</commentary>
  </example>
model: inherit
color: orange
tools: ["Read", "Write", "Edit", "Bash", "Glob", "Grep", "WebFetch", "LSP", "TaskCreate", "TaskUpdate", "TaskList"]
---
```

The agent body is a complete domain expert prompt — it does NOT have a `skills` field in the frontmatter. Instead, it references skill files MANUALLY in the body:

```markdown
## Rust Knowledge Skill Reference Files

This agent has access to comprehensive reference documentation via the rust-knowledge skill.
Use the Read tool to access these files when needed:

| Topic | Path |
|-------|------|
| Ownership & Borrowing | `$CLAUDE_PLUGIN_ROOT/skills/rust-knowledge/references/ownership.md` |
| Error Handling | `$CLAUDE_PLUGIN_ROOT/skills/rust-knowledge/references/error-handling.md` |
```

The agent reads these files using the Read tool ON DEMAND when it needs them — they are NOT auto-loaded.

### Agent + Skill Relationship

No auto-loading. The agent frontmatter does NOT have a `skills` field. The relationship is:

- `SKILL.md` = the "index" file — a quick reference that also lists all sub-references and scripts
- `references/*.md` = deeper domain knowledge files
- `scripts/*.sh` = runnable utilities

The agent body contains a table of WHEN to read each reference file:

```markdown
**When to read reference files:**
- Before implementing error handling → read `error-handling.md`
- Before writing async code → read `async-rust.md`
- Before writing tests → read `testing.md`
- When troubleshooting build issues → read `troubleshooting.md`
```

The SKILL.md itself is also loaded on demand if the agent needs a quick overview. The agent decides which references to read based on the task at hand.

### User-Invocable Skills (Entry Points)

The command file `commands/rust-developer.md` has:

```yaml
---
description: Comprehensive Rust development assistance...
argument-hint: Task description
allowed-tools: ["Task"]
---
```

Its body contains a SINGLE instruction: immediately invoke the `rust-developer` agent via the Task tool:

```markdown
**IMMEDIATELY invoke the agent using the Task tool:**
Task tool parameters:
- subagent_type: "rust-developer:rust-developer"
- description: "Rust development assistance"
- prompt: <the user's request from $ARGUMENTS>

Do not process the request yourself. Launch the agent and let it handle the task.
```

The command is a pure pass-through with zero logic. Its sole job is to translate the user invocation into a Task call with the correct `subagent_type`.

The command does NOT have an `agent:` field in its frontmatter (unlike get-shit-done's plan-phase command). Instead, the body text explicitly instructs Claude to use the Task tool.

### The Handoff Pattern (traced)

Trace for `/rust-developer "create a CLI tool"`:

1. **User types** `/rust-developer "create a CLI tool"`
2. **Claude loads** `commands/rust-developer.md`. The command body says: immediately call `Task(subagent_type="rust-developer:rust-developer", prompt="create a CLI tool")`.
3. **Claude (orchestrator) spawns** the `rust-developer` agent via Task tool.
4. **`rust-developer` agent** receives the task. It reads the agent body (which is its full operating manual). It decides which reference files to read based on the request:
   - Will reference `$CLAUDE_PLUGIN_ROOT/skills/rust-knowledge/references/cargo.md` for Cargo setup
   - Will `WebFetch` official Rust Book docs to verify patterns before implementing
5. **Agent works autonomously**, produces files, runs cargo commands, and returns a structured summary to the orchestrator.
6. **Orchestrator** receives the result and presents it to the user.

There is NO intermediate workflow file. The command → agent transition is direct. All decision-making logic lives inside the agent body.

### Skill Categories Found

Two distinct categories:

**1. Command** (`commands/{plugin-name}.md`) — User-invocable entry point
- Thin pass-through
- Only allowed tool: `Task`
- Body: single instruction to spawn the agent
- No logic, no routing, no state

**2. Agent** (`agents/{plugin-name}.md`) — The actual implementer
- Deep domain knowledge embedded in body
- Explicit instructions for WHEN to read which reference files
- Full operating protocol including code review process, output format, best practices
- References skill files via `$CLAUDE_PLUGIN_ROOT` path variable

**3. Skill** (`skills/{skill-name}/SKILL.md`) — On-demand reference index
- Not auto-loaded
- The agent reads this when it needs an overview or quick reference
- Indexes deeper `references/*.md` files

**4. Reference files** (`skills/{skill-name}/references/*.md`) — Deep domain knowledge
- Read by agent ON DEMAND based on task requirements
- Narrow focus per file (one topic per file)

**5. Scripts** (`skills/{skill-name}/scripts/*.sh`) — Executable utilities
- Agent runs these via Bash tool
- Pre-built, validated patterns (create project, run tests, format code)

**6. Assets** (`skills/{skill-name}/assets/`) — Template files
- Cargo.toml templates, config files, etc.
- Agent copies these when scaffolding new projects

### Output Contracts

Defined inside the agent body as a markdown section:

```markdown
## Output Format

When completing tasks, provide:
1. **Summary**: What was accomplished
2. **Files Changed**: List of files created or modified
3. **Commands Run**: Any cargo commands executed
4. **Verification**: Confirmation that Rust docs were consulted
5. **Next Steps**: Recommendations for follow-up actions
```

No machine-readable format — purely human-readable markdown. The orchestrator (which is just Claude after spawning) does not parse structured tokens; it receives the prose and presents it to the user.

### Routing Mechanism

There is no routing. The pattern is:

- User command → spawn agent
- Agent runs → returns to user

No multi-agent orchestration. No conditional spawning based on flags. No revision loops. No workflow files. The command is a one-shot launcher; the agent handles everything from there.

---

## Synthesis: The Correct Pattern

### How Many Agent Types and What Are Their Scopes

**get-shit-done** uses a **multi-agent swarm** pattern with 11 specialized agents, each scoped to a single operation in a pipeline:

- Research agents (domain research, codebase mapping)
- Planning agents (roadmapper, planner, plan-checker)
- Execution agents (executor, verifier, integration-checker)
- Support agents (debugger, synthesizer)

**claude-code-plugins** uses a **single expert agent** pattern: one agent per domain (Rust, Swift, SQL, PowerShell), scoped to ALL operations within that domain.

For the-rewrite-room plugin (which involves planning, research, execution within a specific domain), the correct model depends on scope:

- If the plugin needs a PIPELINE of operations (research → plan → execute → verify), use the multi-agent pipeline pattern from get-shit-done.
- If the plugin is a domain expert that handles user requests end-to-end within that domain, use the single-agent pattern from claude-code-plugins.

### How Agent Skills (Auto-Loaded) Relate to Command Skills (User-Invocable)

**Neither repo uses auto-loading via a `skills` field.**

In get-shit-done: agents reference workflow files and reference files using `@` inline references in their body. These load when the agent context is initialized.

In claude-code-plugins: agents reference skill files explicitly in a lookup table in their body with conditional loading instructions ("Before implementing X, read reference Y").

**The pattern in both repos**: Skills are NOT auto-loaded into agent context. They are read on-demand by the agent using the Read tool when the task requires that domain knowledge.

**For the-rewrite-room**: Do NOT use a `skills` field in agent frontmatter expecting auto-loading. Instead, include a reference table in the agent body with conditional load instructions.

### What Goes in an Agent's Auto-Loaded Skill vs a Command Skill

There is no "auto-loaded skill" in either repo. The correct framing:

**Command file** contains:
- Invocation metadata (name, description, argument-hint, allowed-tools)
- How to pass `$ARGUMENTS` to the agent
- Whether to use direct `agent:` field or Task tool call in body

**Agent file** contains:
- Full operating protocol (role, responsibilities, step-by-step execution flow)
- Reference lookup table (which skill files to read for which sub-tasks)
- Output format contracts
- Success criteria

**SKILL.md** (the skill index file) contains:
- Quick reference (essential commands, common patterns)
- Index of deeper reference files (what each covers)
- Links to scripts and assets

**`references/*.md`** files contain:
- Deep domain knowledge for one specific topic
- Code examples verified against official documentation
- Anti-patterns to avoid

### How the Decision Tree / Routing Is Encoded

**get-shit-done**: Routing lives in WORKFLOW files (`get-shit-done/workflows/*.md`). Commands are thin; they load the workflow via `@`. The workflow is a numbered step-by-step process with conditional logic that tells the orchestrator (Claude) which agent to spawn under which conditions.

```text
Command → loads Workflow → Workflow tells orchestrator when/which agents to spawn → Agents return structured tokens → Workflow routes based on return tokens
```

**claude-code-plugins**: No routing. Decision tree does not exist. Command spawns agent, agent handles everything.

```text
Command → spawns Agent → Agent handles task end-to-end
```

**For the-rewrite-room**: If the plugin needs conditional orchestration (e.g., research on first run, skip if cached; verify after planning; revision loop on failure), use a workflow file loaded by the command via `@`. If the plugin is a single expert handling user requests, embed all logic in the agent body.

### The Exact Invocation Flow from User Command to Agent Execution

**Simple pattern** (claude-code-plugins):

```text
1. User: /rust-developer "create CLI"
2. Claude loads: commands/rust-developer.md
3. Command body: Task(subagent_type="rust-developer:rust-developer", prompt="create CLI")
4. Agent (rust-developer.md) executes:
   a. Reads task from prompt
   b. Identifies which reference files to load (Read tool, on-demand)
   c. WebFetch official docs to verify patterns
   d. Implements, runs cargo, returns structured output
5. User sees result
```

**Pipeline pattern** (get-shit-done):

```text
1. User: /gsd:plan-phase 3
2. Claude loads: commands/gsd/plan-phase.md
3. Command loads via @: workflows/plan-phase.md (the decision tree)
4. Orchestrator follows plan-phase.md workflow:
   a. Initialize: call gsd-tools.cjs to load state
   b. Check flags, existing artifacts
   c. IF research needed: Task(subagent_type="gsd-phase-researcher", prompt=...)
      Agent returns: ## RESEARCH COMPLETE
   d. Task(subagent_type="gsd-planner", prompt=...)
      Agent returns: ## PLANNING COMPLETE
   e. IF checker enabled: Task(subagent_type="gsd-plan-checker", prompt=...)
      Agent returns: ## VERIFICATION PASSED or ## ISSUES FOUND
   f. IF issues: Task(subagent_type="gsd-planner", prompt=... with revision_context)
      Revision loop max 3 iterations
5. User sees planning result with wave structure
```

### Key Implementation Rules for the-rewrite-room

1. **Command file is always thin**: Its job is to load context and trigger the agent. No business logic in commands.

2. **Agent body IS the operating manual**: Every step the agent takes, every deviation rule, every output format lives in the agent `.md` file body.

3. **Workflow files carry orchestration logic**: If the plugin needs multi-agent coordination or conditional routing, extract that into a workflow file and load it from the command via `@~/.claude/{plugin}/workflows/{name}.md`.

4. **Skill files are read on-demand**: Agents include a lookup table (Markdown table) in their body that maps task types to reference files. The agent reads the relevant files using the Read tool when needed, not upfront.

5. **Output contracts are structured tokens**: Agents return machine-parseable headings (`## PLAN COMPLETE`, `## CHECKPOINT REACHED`) when the orchestrator needs to branch on the result. Use prose output only when the result goes directly to the user.

6. **`$CLAUDE_PLUGIN_ROOT`** is the runtime variable pointing to the installed plugin directory. Use it in agent bodies to construct absolute paths to reference files and scripts.

7. **`agent:` frontmatter field in commands**: This field causes Claude to hand off to the named agent directly (without a Task tool call in the body). It is an alternative to having the body explicitly call `Task(subagent_type=...)`. Use either — both work.

8. **Subagent type syntax**: `"plugin-name:agent-name"` — the plugin name prefix is required when spawning agents from a plugin context. Example: `rust-developer:rust-developer`.
