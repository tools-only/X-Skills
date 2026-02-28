---
name: create-agent
description: >
  This skill should be used when the user asks to "create an agent", "make an agent",
  "write an agent", "build a subagent", "add an agent to a plugin", "design an autonomous agent",
  "generate an agent file", "write a system prompt for an agent", "what frontmatter does an agent need",
  "create a specialized agent", or needs help with agent triggering conditions, expert personas,
  tool restrictions, or the <example> block description format for agents.
argument-hint: "[name] - or leave empty to interview"
---

# Agent Creator

Generate well-structured Claude Code agents — markdown files with YAML frontmatter that delegate
complex multi-step work to autonomous subprocesses with isolated context windows.

**Agents vs Skills — know the difference before generating:**

- **Agents** run in isolated context, have second-person system prompts ("You are..."), use `<example>` XML blocks in descriptions, and are spawned via the Task tool
- **Skills** inject inline into the current conversation, use imperative body instructions for Claude to follow, and route via description matching on trigger phrases

## Phase 0: Fetch Current Documentation

**Before generating**, retrieve the latest agent documentation:

```
Use Task tool with subagent_type=claude-code-guide:
"List all current frontmatter options for Claude Code agents (AGENT.md), including
tools, disallowedTools, model, color, permissionMode, isolation, background, maxTurns,
skills, memory, and any new fields added recently."
```

Capture any fields not in `${CLAUDE_PLUGIN_ROOT}/skills/create-agent/references/agent-frontmatter.md`. If nothing new, proceed. If the subagent fails, proceed — current references are sufficient.

Proceed to Phase 1.

## Phase 1: Understand Requirements

Parse `$ARGUMENTS` for hints. Gather:

1. **Domain & purpose** — What problem does this agent solve?
2. **Expert persona** — What specialist identity should it embody?
3. **Trigger conditions** — When should Claude delegate to this agent? What user messages activate it?
4. **Proactive vs reactive** — Should it fire automatically after events (e.g., after code is written), or only on explicit request?
5. **Tool access** — What tools are actually needed? Least-privilege: an analysis agent doesn't need Write.
6. **Context isolation** — Does it generate heavy output? Should it run in a worktree (`isolation: worktree`)?

Proceed to Phase 2 once domain, trigger conditions, and proactive intent are established.

## Phase 2: Generate

Apply throughout: second-person for system prompt body, intensional over extensional reasoning,
minimum viable frontmatter.

### Step 1 — Choose identifier

Naming rules (enforced by `validate_agent.py`):

- 3–50 characters, lowercase letters/numbers/hyphens only
- Must start and end with alphanumeric
- Avoid generic terms: `helper`, `assistant`, `agent`

Good: `code-reviewer`, `test-generator`, `api-docs-writer`
Bad: `ag` (too short), `-start` (leading hyphen), `my_agent` (underscore)

### Step 2 — Write frontmatter

Read `${CLAUDE_PLUGIN_ROOT}/skills/create-agent/references/agent-frontmatter.md` for the full
field catalog, color semantics, model options, tool selection framework, and execution modifiers.

**Required:** `name`, `description`
**Always set:** `model: inherit` (unless specific model capability needed), `color` (visual ID in UI)

**Intensional rule for `tools`:** Restrict to the minimum needed because agents run autonomously —
over-permission has no human in the loop to catch it. `["Read", "Grep", "Glob"]` for analysis.
Add `Write` for generation. Add `Bash` only when shell execution is essential, never by default.

### Step 3 — Write description with `<example>` blocks

Agent descriptions use a unique two-part format: a trigger statement + XML example blocks.
The routing model uses this structure to decide when to spawn the agent.

**Format:**

```
Use this agent when [trigger conditions]. Examples:

<example>
Context: [Situation description]
user: "[User request]"
assistant: "[How Claude responds before delegating to agent]"
<commentary>
[Why this agent should trigger here]
</commentary>
</example>
```

**Intensional rules for the description field:**

- 2–4 `<example>` blocks — single-example descriptions miss synonym trigger coverage
- Cover different phrasings of the same intent — routing matches on token patterns across examples
- Include proactive examples when the agent should fire after events, not just on explicit request
- `<commentary>` must explain routing reasoning, not just restate the user message
- Specify when NOT to use if ambiguity with other agents exists — prevents mis-routing

For proactive agents, the description needs a two-turn assistant pattern showing an event followed
by delegation. See `references/agent-frontmatter.md` for the exact format.

### Step 4 — Write system prompt

The markdown body (after `---`) becomes the agent's system prompt. Write entirely in second person,
addressing the agent directly. This is the critical authoring difference from skills: agents need
a persona and process, not instructions for Claude to follow.

**Standard structure:**

```markdown
You are [role] specializing in [domain].

**Your Core Responsibilities:**
1. [Primary responsibility]
2. [Secondary responsibility]

**Process:**
1. [Step — imperative]
2. [Next step]

**Quality Standards:**
- [Standard]

**Output Format:**
[Structure and content of what to return]

**Edge Cases:**
- [Situation]: [How to handle]
```

**Intensional rules for the system prompt:**

- Persona first — the expert identity shapes all downstream decisions; establish it in the first sentence
- Process steps prevent "winging it" on complex tasks; each step is an explicit decision boundary
- Output format is non-negotiable — callers need predictable structure to consume results
- Define edge cases in the system prompt; discovered-at-runtime errors cost retries

Keep under 3,000 words. Detailed domain reference belongs in `references/` preloaded via the
`skills:` frontmatter field, not embedded directly in the system prompt.

See `${CLAUDE_PLUGIN_ROOT}/skills/create-agent/examples/proactive-code-reviewer.md` for a
complete working example demonstrating the proactive trigger pattern.

### Step 5 — Script opportunity scan

Read `${CLAUDE_PLUGIN_ROOT}/skills/create-skill/references/script-patterns.md` and apply the
five signal patterns to every step in the agent's system prompt:

| Signal | Question | If yes → |
|--------|----------|----------|
| **Repeated Generation** | Does any step produce the same structure across invocations? | Parameterized script in `scripts/` |
| **Unclear Tool Choice** | Does any step combine tools in a fragile sequence? | Script the procedure |
| **Rigid Contract** | Can you write `--help` text for this step right now? | CLI candidate |
| **Dual-Use Potential** | Would a user run this step from the terminal independently? | Design as proper CLI |
| **Consistency Critical** | Must this step produce identical output for identical inputs? | Script — never LLM generation |

### Step 6 — Check delegation

Scan existing agents and skills before finalizing:

- Does an existing agent cover this domain? Extend it, or tighten scope of the new one
- Are there skills or reference files to preload via `skills:` frontmatter for domain knowledge?
- Are there commands or MCPs this agent should delegate sub-tasks to?

**Always use fully qualified names:**

- `Task: subagent_type=plugin-dev:agent-creator` (not just "agent-creator")
- `Skill: claude-skills:create-skill` (not just "create-skill")

### Step 7 — Validate

When creating a new agent file:

```bash
python3 ${CLAUDE_PLUGIN_ROOT}/skills/create-agent/scripts/validate_agent.py <agent-file> --output json
```

Exit 0 = proceed to Phase 3. Exit 1 = parse the `errors` array; each entry has `field`, `message`,
`severity`. Resolve all `critical` and `major` items before writing to disk.

## Phase 3: Deliver

### Output Paths

| Scope | Location |
|-------|----------|
| User agent (global) | `~/.claude/agents/<name>.md` |
| Project agent | `.claude/agents/<name>.md` |
| Plugin agent | `<plugin-root>/agents/<name>.md` |

Agents in `agents/` are auto-discovered — no registration needed. Plugin agents are namespaced
automatically as `plugin-name:agent-name`.

### Initialize agent file (optional scaffold)

When creating from scratch:

```bash
python3 ${CLAUDE_PLUGIN_ROOT}/skills/create-agent/scripts/init_agent.py <name> --path <agents-dir>
```

Exit 0 = file created with placeholders, proceed to fill content.
Exit 1 = naming collision; ask user to rename or confirm overwrite.

### Explain Your Choices

Present the generated agent with brief rationale:

- **What you set and why** — "Set `model: sonnet` because this agent performs complex multi-file reasoning"
- **What you excluded and why** — "Left `isolation` unset; no git state management needed"
- **Tools selected and why** — explicitly justify each tool; undefended tool access is a design smell

### Write and Confirm

Before writing:

```
Writing to: [path]
This will [create new / overwrite existing] file.
Proceed?
```

### After Creation

Summarize:

- Name and file path
- When it triggers (key trigger conditions)
- Tools granted and why
- Suggested test scenario

Proceed to Phase 4.

## Phase 4: Evaluate

| Dimension | Criteria |
|-----------|----------|
| **Clarity (0-10)** | System prompt unambiguous, objective and persona clear |
| **Trigger Precision (0-10)** | Description + examples cover intended trigger space, not broader |
| **Efficiency (0-10)** | System prompt token economy — maximum guidance per token |
| **Completeness (0-10)** | Covers domain requirements; output format defined; edge cases addressed |
| **Safety (0-10)** | Tools restricted to minimum needed; no runaway permission grants |

**Target: 9.0/10.0.** If below, refine once addressing the weakest dimension, then deliver.

## Validation Checklist

**Structure:**
- [ ] File is `<name>.md` in an `agents/` directory
- [ ] Valid YAML frontmatter with `name` and `description`
- [ ] Markdown body is present and substantial

**Description Quality:**
- [ ] Starts with "Use this agent when..."
- [ ] Contains 2–4 `<example>` blocks
- [ ] Each example has Context, user, assistant, and `<commentary>`
- [ ] Covers different trigger phrasings (synonym coverage)
- [ ] Proactive example included if agent should fire after events

**System Prompt Quality:**
- [ ] Written in second person ("You are...", "You will...")
- [ ] Has clear persona/role statement as first sentence
- [ ] Process steps are numbered and imperative
- [ ] Output format is defined
- [ ] Edge cases addressed
- [ ] Under 3,000 words; domain detail offloaded to references/

**Frontmatter:**
- [ ] `model: inherit` unless specific model needed
- [ ] `color` set and semantically meaningful
- [ ] `tools` restricted to minimum needed
- [ ] No TODO placeholders remaining

## Error Handling

| Issue | Action |
|-------|--------|
| Unclear domain | Ask: what does success look like for this agent? |
| Scope too broad | Split into 2–3 focused agents with non-overlapping trigger conditions |
| Conflicts with existing agent | Note overlap; narrow triggering scope or extend the existing one |
| Naming collision | Ask to rename or confirm overwrite |
| System prompt too long | Move domain detail to `references/` and preload via `skills:` field |
| Vague trigger conditions | Ask for 3 concrete user messages that should activate this agent |
