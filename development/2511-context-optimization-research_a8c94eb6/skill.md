# Context Window Optimization Research

Research conducted: 2026-02-24
Sources: docs.anthropic.com (official Claude Code documentation)

---

## 1. Current Context Loading Behavior

### CLAUDE.md and rules files

- CLAUDE.md files in the directory hierarchy **above** the working directory are loaded **in full at launch**.
- CLAUDE.md files in **child directories** (subtrees below cwd) load **on demand** — only when Claude reads files in those subdirectories.
- `.claude/rules/*.md` files without a `paths` frontmatter field are loaded **unconditionally** at launch.
- `.claude/rules/*.md` files **with** a `paths` frontmatter field are **conditional** — loaded only when Claude works with files matching the glob pattern.
- Auto memory (`~/.claude/projects/<project>/memory/MEMORY.md`): only the **first 200 lines** load into the system prompt at session start. Topic files (`debugging.md`, `patterns.md`, etc.) are **not loaded at startup** — Claude reads them on demand using standard file tools.

### Skills context loading

Skill descriptions are loaded into context at session start so Claude knows what skills are available. Full skill content loads only when the skill is invoked.

- There is a character budget for skill descriptions: **2% of context window**, with a fallback of **16,000 characters**. Skills exceeding the budget are excluded. Overridable with `SLASH_COMMAND_TOOL_CHAR_BUDGET` env var.
- Skills with `disable-model-invocation: true` are excluded from context entirely — their description is not loaded. Only loadable by explicit user `/skill-name` invocation.
- Skills with `user-invocable: false` have their description in context (Claude-invocable), but are hidden from the `/` menu.

| Frontmatter                      | Description in context | Full content loaded         |
| :------------------------------- | :--------------------- | :-------------------------- |
| (default)                        | Yes (at session start) | On invocation               |
| `disable-model-invocation: true` | No                     | Only on user `/name` invoke |
| `user-invocable: false`          | Yes (at session start) | On invocation               |

### Subagents context loading

- Subagents are loaded at session start (file scanning). If you add a subagent file mid-session, restart or use `/agents` to load it immediately.
- Each subagent runs in its **own context window** — it does not inherit the parent conversation context.
- Subagents receive only their own system prompt (plus basic environment details like working directory), **not** the full Claude Code system prompt.
- Subagents do **not** inherit skills from the parent conversation. Skills must be listed explicitly in the `skills` frontmatter field.

---

## 2. `skills` Field in Agent Frontmatter

**Confirmed: exists.** Documented in the `sub-agents.md` page under "Supported frontmatter fields".

```yaml
skills:
  - api-conventions
  - error-handling-patterns
```

Behavior (from official docs):
- The **full content** of each skill is injected into the subagent's context at startup — not just made available for invocation.
- Subagents do not inherit skills from the parent conversation; skills must be listed explicitly.
- This is the inverse of `context: fork` in a skill: with `skills` in a subagent, the subagent controls the system prompt and loads skill content. With `context: fork` in a skill, the skill content is injected into the specified agent.

This means `skills` in agent frontmatter is **always-loaded** (at subagent startup) — there is no lazy loading within a subagent's own context via this field. It is a mechanism to pre-inject knowledge, not to conditionally load it.

---

## 3. `hooks` Field in Agent Frontmatter

**Confirmed: exists.** Documented in the `sub-agents.md` page under "Supported frontmatter fields" and in the "Define hooks for subagents" section.

```yaml
hooks:
  PreToolUse:
    - matcher: "Bash"
      hooks:
        - type: command
          command: "./scripts/validate-readonly-query.sh"
  PostToolUse:
    - matcher: "Edit|Write"
      hooks:
        - type: command
          command: "./scripts/run-linter.sh"
```

Behavior (from official docs):
- Hooks defined in agent frontmatter run **only while that specific subagent is active** and are cleaned up when it finishes.
- All hook events are supported: `PreToolUse`, `PostToolUse`, `Stop` (converted to `SubagentStop` at runtime).
- Hooks can block tool execution: a script exiting with code 2 blocks the operation and feeds stderr back to Claude.
- `SessionStart` is NOT listed as a supported hook event for agent frontmatter. The supported events documented for agent frontmatter are `PreToolUse`, `PostToolUse`, and `Stop`.

There is also a **project-level hooks** mechanism in `settings.json` that responds to subagent lifecycle events:
- `SubagentStart` — fires when a subagent begins execution (matcher: agent type name)
- `SubagentStop` — fires when a subagent completes (matcher: agent type name)

These project-level hooks can be used to run setup/teardown scripts around subagent execution from the main session, but they do not load context into the subagent itself.

---

## 4. Thin-Wrapper-Agent Pattern

### Evidence found

The documentation does not explicitly name a "thin-wrapper-agent" or "topic-routing" pattern. However, the following primitives exist that could compose such a pattern:

**Mechanism 1: `disable-model-invocation: true` on skills (lazy loading)**

Skills with `disable-model-invocation: true` are excluded from context entirely until explicitly invoked. This is the primary documented mechanism for on-demand skill loading. A lightweight "router" Claude.md could avoid pre-loading heavy skills by using this flag, then invoking them via `/skill-name` when appropriate.

**Mechanism 2: Subagents with explicit `skills` field**

A subagent can preload a specific set of skills relevant to its domain. A thin orchestrator agent could delegate to domain-specific subagents, each of which loads only its own skill set. This approximates topic-routing: route task to domain subagent → subagent loads domain skills → subagent executes → returns summary.

**Mechanism 3: `context: fork` in skills + `agent` field**

A skill can run in an isolated subagent context via `context: fork`. Combined with `agent: Explore` or a named custom agent, this allows a lightweight skill to delegate into a specialized context. The skill content becomes the subagent's task; the agent field determines what tools and model it uses.

**Mechanism 4: Conditional rules via `paths` frontmatter**

`.claude/rules/*.md` files with `paths` frontmatter only load when Claude works with matching files. This is a built-in lazy-loading mechanism for rule context.

### Pattern that does NOT exist (not found)

- No documented `SessionStart` hook for agent frontmatter that dynamically loads skills based on the user's first message. The `SubagentStart` hook exists in `settings.json` but fires in the **main session**, not inside the subagent — it cannot inject skills into the subagent's context.
- No documented mechanism where a thin agent inspects a request and programmatically selects which skills to load before responding. Skill loading at subagent startup is declarative (listed in `skills` field), not dynamic.
- No "routing orchestrator" pattern described that uses hooks to conditionally load skill sets. The hooks documented for agent frontmatter are tool validation hooks (`PreToolUse`, `PostToolUse`), not context-injection hooks.

---

## 5. Concrete Optimization Opportunities Identified

### A. Use `disable-model-invocation: true` on heavy skills

Skills that are rarely needed but large (e.g., detailed reference docs, API specs) should use `disable-model-invocation: true`. This removes their description from the context budget entirely. The user invokes them explicitly when needed. This is the primary documented mechanism for lazy-loading.

### B. Decompose large skills into SKILL.md + supporting files

The official docs recommend keeping `SKILL.md` under 500 lines. Move detailed reference material into separate files (`reference.md`, `examples.md`) referenced from `SKILL.md`. Claude only reads those files on demand when needed — they do not load at skill invocation.

### C. Use `.claude/rules/` with `paths` frontmatter

Instead of loading all rules unconditionally, scope rules to file patterns. Rules without `paths` load for every session. Rules with `paths` load only when Claude works with matching files. This is built-in conditional loading.

### D. Use auto memory topic files instead of monolithic CLAUDE.md

Auto memory loads only the first 200 lines of `MEMORY.md`. Detailed content in topic files (`debugging.md`, `api-conventions.md`) is NOT loaded at startup — Claude reads them on demand. This mirrors the skill supporting-files pattern.

### E. Domain subagents with scoped skill sets

Instead of loading all skills into the main conversation, create subagents per domain (e.g., `api-developer`, `test-writer`), each with only the skills relevant to that domain in their `skills` field. The main conversation stays lean; domain context loads only when the subagent is delegated a task.

### F. `context: fork` for isolated, high-output tasks

Tasks that produce large output (test runs, documentation fetches, log analysis) should use `context: fork` or explicit subagent delegation. The verbose output stays in the subagent's context; only the summary returns to the main conversation.

### G. CLAUDE.md hierarchy loading is already lazy for subdirectories

Child-directory CLAUDE.md files load on demand. Structure large repos so detailed context lives in subdirectory CLAUDE.md files, not in the root CLAUDE.md. The root file stays lean.

---

## 6. Sources Consulted

| URL | Access date | Status |
| :-- | :---------- | :----- |
| `https://docs.anthropic.com/en/docs/claude-code/settings.md` | 2026-02-24 | Fetched (full content, 110KB persisted) |
| `https://docs.anthropic.com/en/docs/claude-code/memory.md` | 2026-02-24 | Fetched (full content returned inline) |
| `https://docs.anthropic.com/en/docs/claude-code/sub-agents.md` | 2026-02-24 | Fetched (full content returned inline) |
| `https://docs.anthropic.com/en/docs/claude-code/skills.md` | 2026-02-24 | Fetched (full content returned inline) |
| `mcp__claude_ai_Ref__ref_search_documentation` | 2026-02-24 | Search returned errors (error_hash 83878a32e6b1) — searches not usable |

The settings.md content was too large to return inline (110KB) and was saved to a persisted tool-result file. Key facts from that page were not directly read due to size. Settings page facts cited in this document are from memory.md and sub-agents.md cross-references, or from the skills.md page (which covers `SLASH_COMMAND_TOOL_CHAR_BUDGET` and the 2% / 16,000-character budget).

**Note on settings.md**: The `SLASH_COMMAND_TOOL_CHAR_BUDGET` env var and the 2% context-window / 16,000-character fallback budget for skill descriptions came from the skills.md page's Troubleshooting section, not from settings.md directly.
