# Frontmatter Options & Patterns Reference

Load before writing frontmatter in Phase 2, Step 2. Contains the full field catalog,
description patterns, execution modifiers, tool selection framework, and progressive
disclosure patterns.

---

## Common Frontmatter Options

```yaml
---
name: identifier                    # Required for skills
description: >                      # How it's described/triggered
  [See description patterns below]

# Execution modifiers
model: sonnet                       # haiku (fast), sonnet (balanced), opus (complex)
context: fork                       # Run in isolated sub-agent, preserves main context
agent: Explore                      # Route to specialized agent (Explore, Plan, custom)

# Tool access
allowed-tools:                      # Restrict available tools
  - Read
  - Grep
  - Bash(git:*)

# Lifecycle hooks (optional)
hooks:
  PreToolUse:
    - command: "validation-script.sh"
  PostToolUse:
    - command: "cleanup.sh"

# Behavior modifiers
user-invocable: true                # Show in /command menu (default true)
disable-model-invocation: true      # Prevent programmatic invocation (commands only)
argument-hint: [arg1] [arg2]        # Document expected arguments
---
```

---

## Description Patterns

### For Skills (auto-triggered) — principles

- **Third-person framing is a routing signal, not a stylistic choice.** The routing model evaluates the description as a triggering condition. First-person ("Use this skill when...") reads as an instruction to execute. Third-person ("This skill should be used when...") reads as a condition to test. The framing changes how the model interprets the field.
- **Quoted phrases must be verbatim user speech.** Routing matches on literal token patterns. Write the exact words a user would type, not paraphrases: `"create a hook"` triggers correctly; `"hook creation workflows"` may not.
- **The description is always in context, even when the skill isn't active.** Every session pays the token cost of every skill's description. Density matters: cover more trigger patterns in fewer words. Avoid restating the skill name or explaining what skills are.
- **Cover the naive phrasing.** A user who doesn't know this skill exists won't search for it by name — they'll describe their problem in plain language. Include the phrasing someone would use who has never heard of this skill.
- **3–5 trigger phrases minimum.** Single-phrase descriptions have high miss rates. Varied phrases improve routing coverage across synonym space.
- **Use `>` scalar, not `|`.** Folded scalar (`>`) collapses newlines to spaces, producing a single continuous string — correct for descriptions. Literal scalar (`|`) preserves newlines, which can create unexpected whitespace when parsed.

```yaml
# Correct — third-person, verbatim phrases, folded scalar
description: >
  This skill should be used when the user asks to "create a hook",
  "add validation", "implement lifecycle automation", or mentions
  pre/post tool events.

# Wrong — vague, no trigger phrases, not third-person
description: Provides guidance for hooks.
```

### For Commands (user-invoked) — principles

- **Verb-first, under 60 chars.** The description appears as a single scannable line in the `/` menu — treat it as a menu label, not a sentence.
- **Describe the action, not the tool.** "Fix GitHub issue by number" orients by outcome. "GitHub issue fixer" orients by tool name. Users scan for what they want to accomplish.

```yaml
description: Fix GitHub issue by number
description: Review code for security issues
description: Deploy to staging environment
```

---

## Execution Modifiers

Use these when the default behavior isn't sufficient:

- **`context: fork`** — Run in isolated sub-agent. Use for heavy workflows that would pollute main context, or when you need clean separation. Requires `agent:`.

- **`agent: [type]`** — Route to a specialized agent. Examples: `Explore` for codebase search, `Plan` for architecture decisions, or custom agents you've defined.

- **`model: [level]`** — Override the model. Valid values: `haiku` (fast, cheap, simple tasks), `sonnet` (balanced default), `opus` (complex reasoning). Omit to inherit from the current conversation.

- **`hooks`** — Run scripts before/after tool use, scoped to this skill's lifecycle. Useful for validation, logging, or side effects.

- **`disable-model-invocation: true`** — Prevent Claude from auto-loading this skill. Use for skills you want to invoke manually only (commands only).

- **`user-invocable: false`** — Hide from the `/` command menu. Use for background-knowledge skills that should trigger automatically but not appear as slash commands.

---

## Tool Selection

Default generous, restrict only when needed. The principle: restrict tools that have destructive or side-effect potential, not tools that are read-only or purely generative.

| Tier | Tools | Why |
|------|-------|-----|
| **Always allow** | Read, Grep, Glob | Read-only, no side effects |
| **Usually allow** | Edit, Write, WebSearch, WebFetch, Task | Core work tools; restrict if skill is deliberately read-only |
| **Scope Bash** | `Bash(git:*)`, `Bash(npm:*)`, `Bash(pytest:*)` | Bash is the highest blast-radius tool — scope to known commands |
| **If interactive** | AskUserQuestion | Required any time the skill needs user decisions mid-workflow |
| **If delegating** | Skill | Required to invoke other skills programmatically |
| **If notebooks** | NotebookEdit | Jupyter-specific; omit unless skill touches `.ipynb` files |
| **If plan-gated** | ExitPlanMode, EnterPlanMode | For workflows requiring explicit user approval before execution |

---

## Progressive Disclosure

For complex skills, organize into subdirectories:

```
skill-name/
├── SKILL.md          # Core instructions (keep under 500 lines)
├── scripts/          # Executable code (Python/Bash)
├── references/       # Docs loaded into context as needed
├── examples/         # Working code examples users can copy directly
└── assets/           # Files used in output (templates, icons, fonts)
```

**scripts/** — Deterministic, token-efficient. May be executed without loading into context. Use when the same code is rewritten repeatedly or reliability is critical.

**references/** — Documentation Claude reads while working. Keeps SKILL.md lean. For files >100 lines, include a table of contents. Only load when needed.

**examples/** — Working code examples: complete, runnable scripts, configuration files, template files, real-world usage examples. Users can copy and adapt these directly. Distinct from references (docs) and scripts (utilities).

**assets/** — Files NOT loaded into context. Used in output: templates, images, fonts, boilerplate.

### Pattern 1: High-level guide with references

```markdown
# PDF Processing

## Quick start
Extract text with pdfplumber:
[code example]

## Advanced features
- **Form filling**: See references/forms.md
- **API reference**: See references/api.md
```

Claude loads references only when needed.

### Pattern 2: Domain-specific organization

```
bigquery-skill/
├── SKILL.md (overview and navigation)
└── references/
    ├── finance.md (revenue, billing)
    ├── sales.md (pipeline, opportunities)
    └── product.md (API usage, features)
```

When user asks about sales, Claude only reads sales.md.

### Pattern 3: Variant-based organization

```
cloud-deploy/
├── SKILL.md (workflow + provider selection)
└── references/
    ├── aws.md
    ├── gcp.md
    └── azure.md
```

User chooses AWS → Claude only reads aws.md.
