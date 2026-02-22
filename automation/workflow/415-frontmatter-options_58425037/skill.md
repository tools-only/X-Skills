# Frontmatter Options Reference

Complete reference for auditing skill and command frontmatter. Load this before running
Dimension 2 audits. Every field listed here is the full set of valid options — anything
not on this list is not a valid frontmatter key.

## Fields

### `name` (string)
Lowercase letters, digits, and hyphens only. Max 64 characters. Omit to use the directory
name as the identifier. Prefer short, verb-led names for commands. Namespace by tool when
it aids routing clarity: `gh-address-comments`, `linear-close-issue`.

### `description` (string)
The primary triggering mechanism. Always in context — costs tokens on every session
regardless of whether the skill is active. Use `>` folded scalar (not `|` literal). Must
be third-person for skills ("This skill should be used when..."), verb-first under 60
chars for commands.

### `model` (enum)

| Value | When to use |
|-------|-------------|
| `haiku` | Single-step operations, formatting, simple lookups, anything token-cheap and time-sensitive |
| `sonnet` | Balanced default for most workflows; multi-step reasoning without deep analysis |
| `opus` | Complex multi-step reasoning, architectural decisions, deep analysis across large codebases |
| `inherit` | Default when field is omitted — uses the model of the current conversation |

Audit rule: `opus` is only justified when the task genuinely requires deep reasoning that
`sonnet` would struggle with. Most skills should omit `model:` entirely and inherit.
`haiku` is underused — prefer it for lookup-heavy, low-judgment steps.

### `context` (enum)

| Value | When to use |
|-------|-------------|
| `fork` | Only valid value. Runs the skill in an isolated sub-agent. |

Use `context: fork` when: the skill produces heavy output that would pollute the main
conversation context; the skill needs clean separation from conversation history; the
skill delegates to a specialized agent via `agent:`. Requires `agent:` to be set.

### `agent` (string)

Requires `context: fork`. Routes execution to a specialized agent type.

| Value | Capabilities |
|-------|-------------|
| `Explore` | Codebase search, file discovery, read-only analysis. No write tools. |
| `Plan` | Architecture and design planning. Read-only. Outputs plans for approval. |
| `general-purpose` | Full tool access including edit, write, bash. |
| `<custom-name>` | Any agent defined in `.claude/agents/` in the current project. |

### `allowed-tools` (list)

Restricts which tools the skill can use. Default is unrestricted. Specifying this list
is a security and scope constraint — use it to limit blast radius for sensitive skills.

**Complete tool list:**

| Tool | Category | Notes |
|------|----------|-------|
| `Read` | File ops | Read-only, no side effects |
| `Write` | File ops | Creates/overwrites files |
| `Edit` | File ops | String replacement in existing files |
| `Glob` | Search | Pattern-based file discovery |
| `Grep` | Search | Regex content search |
| `Bash` | Execution | Highest blast radius — scope with patterns |
| `Bash(git:*)` | Execution | Scoped to git commands only |
| `Bash(npm:*)` | Execution | Scoped to npm commands only |
| `Bash(pytest:*)` | Execution | Scoped to pytest commands only |
| `WebFetch` | Web | Fetches a specific URL; distinct from search |
| `WebSearch` | Web | Queries a search engine |
| `Task` | Orchestration | Spawns subagents |
| `AskUserQuestion` | Interaction | Required for any mid-workflow user decision |
| `Skill` | Invocation | Required to invoke other skills programmatically |
| `NotebookEdit` | Notebooks | Jupyter-specific; omit unless skill touches `.ipynb` |
| `EnterPlanMode` | Plan flow | Required for plan-gated workflows |
| `ExitPlanMode` | Plan flow | Required for plan-gated workflows |
| `mcp__<server>__<tool>` | MCP | Any tool from an installed MCP server |

**Audit rule:** `Bash` unrestricted is almost always wrong — scope it. `AskUserQuestion`
must be present if the skill asks the user anything. `Skill` must be present if the skill
invokes another skill by name.

### Tool Selection Framework

The core principle: restrict tools that have destructive or side-effect potential, not
tools that are read-only or purely generative. Over-restriction breaks the skill; under-
restriction is a security and scope risk.

| Tier | Tools | Why | When to restrict |
|------|-------|-----|------------------|
| **Always allow** | `Read`, `Grep`, `Glob` | Read-only, no side effects | Only if skill must be strictly read-only |
| **Usually allow** | `Edit`, `Write`, `WebSearch`, `WebFetch`, `Task` | Core work tools | Restrict if skill is deliberately non-modifying |
| **Scope Bash** | `Bash(git:*)`, `Bash(npm:*)`, `Bash(pytest:*)` | Highest blast radius — scope to known commands | Never allow unrestricted `Bash` unless tool scope is genuinely unknown |
| **Require if interactive** | `AskUserQuestion` | Required any time the skill needs user decisions | Omit only if the skill is fully automated |
| **Require if delegating** | `Skill` | Required to invoke other skills programmatically | Omit if no delegation |
| **Require if notebook** | `NotebookEdit` | Jupyter-specific | Omit unless skill touches `.ipynb` |
| **Require if plan-gated** | `EnterPlanMode`, `ExitPlanMode` | For workflows requiring explicit approval before execution | Omit unless skill has a plan/execute split |

**Gap audit questions:**
- Does the skill need user decisions but lacks `AskUserQuestion`? → add it
- Does the skill invoke another skill but lacks `Skill`? → add it
- Does the skill write files but has no `Edit` or `Write`? → add them
- Does the skill have unrestricted `Bash` when a scoped pattern would work? → scope it
- Does the skill have `Bash(git:*)` but never uses git? → remove it (dead scope)

### `hooks` (object)

Scoped to this skill's lifecycle. Runs scripts before or after specific tool events.

```yaml
hooks:
  PreToolUse:
    - command: "scripts/validate-input.sh"
  PostToolUse:
    - command: "scripts/cleanup.sh"
  Stop:
    - command: "scripts/on-complete.sh"
```

Valid hook events: `PreToolUse`, `PostToolUse`, `Stop`. Use for validation, logging,
side effects, or cleanup that must happen deterministically around tool calls.

### `user-invocable` (boolean)

| Value | Effect |
|-------|--------|
| `true` | Default. Skill appears in the `/` command menu. |
| `false` | Hidden from `/` menu. Skill still triggers automatically via description routing. |

Use `false` for background-knowledge skills that should activate automatically but
shouldn't clutter the command menu.

### `disable-model-invocation` (boolean)

Commands only. Prevents Claude from auto-loading this skill based on its description.
Forces manual invocation only. Default: `false`.

### `argument-hint` (string)

Shown in autocomplete when the user types the command. Documents expected argument
syntax. Examples: `[issue-number]`, `<path-to-skill>`, `[skill|command] [name]`.

Audit rule: any skill that reads `$ARGUMENTS` or `$1`/`$2` should have `argument-hint`
set so users know what to pass.

---

## Dynamic Content Syntax

These substitutions are processed before the skill body reaches Claude.

| Syntax | Resolves to |
|--------|-------------|
| `$ARGUMENTS` | All arguments passed to the skill as a single string |
| `$1`, `$2`, `$3` | Individual positional arguments |
| `@path/to/file` | Contents of the file at that path, loaded inline |
| `@$1` | Contents of the file whose path was passed as the first argument |
| bang + backtick-wrapped command (e.g. `!date`) | Output of executing the command in a shell, injected inline |

Audit rule: skills that accept a file path as input should use `@$1` to load it inline
rather than requiring a separate Read tool call — the injection happens before the model
sees the skill, saving a tool round-trip. The bang-backtick pattern is underused:
real-time data like git branch, file tree, or env vars can be injected without tool calls.
