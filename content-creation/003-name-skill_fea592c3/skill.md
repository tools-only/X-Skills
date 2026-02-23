---
name: create-skill
description: >
  This skill should be used when the user asks to "create a skill", "make a command",
  "generate a prompt", "write a slash command", "build a Claude extension",
  "add a skill to a plugin", "improve skill description", or "write skill frontmatter".
argument-hint: "[skill|command] [name] - or leave empty to interview"
---

# Skill & Command Generator

Generate well-structured skills or slash commands. Both are markdown files with YAML frontmatter—they share the same structure but differ in how they're triggered and described.

## Phase 0: Fetch Current Documentation

**Before generating**, retrieve the latest documentation:

```
Use Task tool with subagent_type=claude-code-guide:
"List all current frontmatter options for skills and commands, including any execution modifiers, model selection, and structural options."
```

Capture any frontmatter fields or options not already listed in `references/frontmatter-options.md`. If nothing new, proceed with current references. If the subagent fails, proceed — current references are sufficient.

## Phase 1: Understand Requirements

Parse `$ARGUMENTS` for type hint. User is often unclear and uninformed on best practices of skill development. Continue to interview and help user using `/brainstorm` skill, if available.

Gather requirements:
1. **Primary objective** — What should this do?
2. **Trigger scenarios** — When should it activate?
3. **Inputs/outputs** — What does it receive and produce?
4. **Complexity** — Simple, standard, or complex?
5. **Execution needs** — Isolated context? Delegated to specialized agent?

Proceed to Phase 2 when at minimum Objective and Trigger Scenarios are established. Remaining dimensions can be resolved during generation.

## Phase 2: Generate

Apply these principles throughout generation: use imperative voice and terse phrasing because every token in a generated skill body costs budget on every invocation, and Claude extrapolates well from precise nudges. Prefer instruction over example — state the rule with its reasoning so it generalizes to every input.

### Step 1 — Choose type

- **Skills:** Trigger-rich, third-person description ("This skill should be used when..."); auto-triggered by routing
- **Commands:** Concise, verb-first description, under 60 chars; user-invoked via `/` menu

### Step 2 — Write frontmatter

Read `${CLAUDE_PLUGIN_ROOT}/skills/create-skill/references/frontmatter-options.md` for the full field catalog, description patterns, tool selection framework, and execution modifiers.

**Intensional over extensional — apply to all generated content.** State the rule directly with its reasoning rather than listing examples that imply the rule. An intensional rule ("quoted phrases must be verbatim user speech *because* routing matches on literal tokens") generalizes to every input the skill will encounter. An extensional approach requires the reader to reverse-engineer the rule — two reasoning hops instead of one, covering only the shape of those specific examples. Since this skill generates instructions that will themselves guide further generation, the quality of reasoning propagates.

### Step 3 — Write body

**Construction rules:**
- State objective explicitly in first sentence
- Use imperative voice ("Analyze", "Generate", "Identify") — no first-person ("I will", "I am")
- Context only when necessary for understanding
- XML tags only for complex structured data
- No "When to Use This Skill" section — body loads only after triggering; routing guidance there is never read by the routing decision
- Avoid headers deeper than H3 — deep nesting signals content that belongs in `references/`, not `SKILL.md`
- Use bang-backtick syntax for dynamic context injection when real-time data (git status, file list, env vars) improves the skill without requiring a tool call

Both skills and commands follow the same body pattern:

```markdown
# Name

Brief overview (1-2 sentences).

## Process
1. Step one (imperative voice)
2. Step two
3. Step three
```

**Dynamic Content:**

| Syntax | Purpose |
|--------|---------|
| `$ARGUMENTS` | All arguments as string |
| `$1`, `$2`, `$3` | Positional arguments |
| `@path/file` | Load file contents |
| `@$1` | Load file from argument |
| Exclamation + backticks | Execute bash command, include output |

### Step 4 — Script opportunity scan

Read `${CLAUDE_PLUGIN_ROOT}/skills/create-skill/references/script-patterns.md` and apply the five signal patterns to every workflow step in the skill being generated:

| Signal | Question | If yes → |
|--------|----------|----------|
| **Repeated Generation** | Does any step produce the same structure with different params across invocations? | Parameterized script in `scripts/` |
| **Unclear Tool Choice** | Does any step combine multiple tools in a fragile sequence naturally expressible as one function? | Script the procedure |
| **Rigid Contract** | Can you write `--help` text for this step right now without ambiguity? | CLI candidate — delegate design to `create-cli` |
| **Dual-Use Potential** | Would a user want to run this step from the terminal, outside the skill workflow? | Design as proper CLI from the start |
| **Consistency Critical** | Must this step produce bit-for-bit identical output for identical inputs? | Script — never LLM generation |

For each identified script candidate:
1. Choose the archetype from `references/script-patterns.md` (init/validate/transform/package/query)
2. If the interface is non-trivial, delegate to `create-cli` skill to design it
3. Scaffold the script in `scripts/` using the Python template from `references/script-patterns.md`
4. Wire it into SKILL.md with: trigger condition, exact invocation, output interpretation

**Wiring rule:** A script reference must state *when* to invoke (trigger condition), *how* to invoke (exact command with flags), and *what to do* with the result (exit code handling, which output fields matter).

### Step 5 — Check delegation

Scan for existing resources before finalizing:

```
Review available: skills, commands, agents, MCPs
For each workflow step, ask: "Do we already have this?"
```

**Common delegation patterns:**
- Git commits → `SlashCommand: /commit`
- Code review → `Skill: /code-review`

**Always use fully qualified names:**
- `Skill: plugin-dev:hook-development` (not just "hook-development")
- `SlashCommand: /plugin-dev:create-plugin` (not just "create-plugin")
- `Task: subagent_type=plugin-dev:agent-creator`

### Step 6 — Validate

When generating a new skill directory (not editing an existing single file):

```bash
python3 ${CLAUDE_PLUGIN_ROOT}/skills/create-skill/scripts/validate_skill.py <skill-directory> --output json
```

Exit 0 = proceed to Phase 3. Exit 1 = parse the `errors` array; each entry has `field`, `message`, `severity`. Resolve all `critical` and `major` items before writing to disk.

### Explain Your Choices

When presenting the generated skill/command to the user, briefly explain:
- **What you set and why** — "Added `context: fork` because this workflow generates heavy output"
- **What you excluded and why** — "Left `model` unset (inherits default), `hooks` omitted (no validation needed)"
- **Add more trigger phrases if routing misses expected inputs**

## Phase 3: Deliver

### Output Paths

| Type | Location |
|------|----------|
| User skill | `~/.claude/skills/<name>/SKILL.md` |
| User command | `~/.claude/commands/<name>.md` |
| Project skill | `.claude/skills/<name>/SKILL.md` |
| Project command | `.claude/commands/<name>.md` |

### Write and Confirm

Before writing:
```
Writing to: [path]
This will [create new / overwrite existing] file.
Proceed?
```

### After Creation

Summarize what was created:
- Name and type
- Path
- How to invoke/trigger
- Suggested test scenario

## Phase 4: Evaluate

Score the generated skill/command:

| Dimension | Criteria |
|-----------|----------|
| **Clarity (0-10)** | Instructions unambiguous, objective clear |
| **Precision (0-10)** | Appropriate specificity without over-constraint |
| **Efficiency (0-10)** | Token economy—maximum value per token |
| **Completeness (0-10)** | Covers requirements without gaps or excess |
| **Usability (0-10)** | Practical, actionable, appropriate for target use |

**Target: 9.0/10.0.** If below, refine once addressing the weakest dimension, then deliver.

## Bundled Scripts

**Initialize a new skill — invoke at Step 1 when creating a new skill directory (not when editing an existing file):**

```bash
python3 ${CLAUDE_PLUGIN_ROOT}/skills/create-skill/scripts/init_skill.py <name> --path <dir> [--resources scripts,references,assets] [--examples]
```

Exit 0 = directory created, proceed to Step 2. Exit 1 = naming collision; ask user whether to overwrite or rename.

**Validate a skill — invoke at Step 6 before delivering:**

```bash
python3 ${CLAUDE_PLUGIN_ROOT}/skills/create-skill/scripts/validate_skill.py <skill-directory> --output json
```

Exit 0 = proceed to Phase 3. Exit 1 = parse `errors` array; resolve all `critical` and `major` items before writing.

**Package for distribution — invoke only when user explicitly requests a distributable file:**

```bash
python3 ${CLAUDE_PLUGIN_ROOT}/skills/create-skill/scripts/package_skill.py <skill-directory> [output-dir]
```

Exit 0 = `.skill` file created at output path. Exit 1 = validation failed; read stdout for details.

## Degrees of Freedom

Match instruction specificity to the task's fragility and variability:

| Level | When to Use | Format |
|-------|-------------|--------|
| **High freedom** | Multiple valid approaches, context-dependent decisions | Text instructions, heuristics |
| **Medium freedom** | Preferred pattern exists, some variation acceptable | Pseudocode, scripts with parameters |
| **Low freedom** | Fragile operations, consistency critical, specific sequence required | Exact scripts, few parameters |

## Quality Standards

**Format Economy:**
- Simple task → direct instruction, no sections
- Moderate task → light organization with headers
- Complex task → full semantic structure

**Balance Flexibility with Precision:**
- Loose enough for creative exploration
- Tight enough to prevent ambiguity

**Remove ruthlessly:** Filler phrases, obvious implications, redundant framing, excessive politeness

## Validation Checklist

Before finalizing a skill or command:

**Structure:**
- [ ] SKILL.md file exists with valid YAML frontmatter
- [ ] Frontmatter has `name` and `description` fields
- [ ] Markdown body is present and substantial
- [ ] Referenced files actually exist

**Description Quality:**
- [ ] Uses third person ("This skill should be used when...")
- [ ] Includes specific trigger phrases users would say
- [ ] Lists concrete scenarios ("create X", "configure Y")
- [ ] Not vague or generic

**Content Quality:**
- [ ] Body uses imperative/infinitive form, not second person
- [ ] Body is focused and lean (1,500–2,000 words ideal, <5k max)
- [ ] Detailed content moved to references/
- [ ] Examples are complete and working
- [ ] Scripts are executable and documented
- [ ] Script opportunities identified via five signal patterns (references/script-patterns.md)
- [ ] Script references in SKILL.md include trigger condition, invocation, output handling
- [ ] Consistency-critical steps are scripted, not left to LLM re-generation

**Progressive Disclosure:**
- [ ] Core concepts in SKILL.md
- [ ] Detailed docs in references/
- [ ] Working code in examples/
- [ ] Utilities in scripts/
- [ ] SKILL.md references these resources

**Testing:**
- [ ] Skill triggers on expected user queries
- [ ] Content is helpful for intended tasks
- [ ] No duplicated information across files
- [ ] References load when needed

## Error Handling

| Issue | Action |
|-------|--------|
| Unclear requirements | Ask clarifying questions |
| Missing context | Request examples or constraints |
| Path issues | Verify directory exists, create with confirmation |
| Type unclear | Default to skill if auto-triggering desired |

---

Execute phases sequentially. Always fetch current documentation first.
