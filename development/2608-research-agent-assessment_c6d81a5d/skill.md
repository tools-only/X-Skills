# Research Agent Assessment — Phase 1 Plugin Creation

**Assessed:** 2026-02-26
**Scope:** All agents in `plugins/plugin-creator/agents/`
**Purpose:** Identify best agent for Phase 1 research tasks in plugin creation workflows

---

## Research Tasks Being Evaluated

1. Search existing plugins/skills for similar functionality (code discovery)
2. Gather domain knowledge about plugin structure, conventions, frontmatter rules
3. Identify architecture patterns from well-structured plugins

---

## Agent Inventory

### 1. plugin-assessor

**File:** `plugins/plugin-creator/agents/plugin-assessor.md`

**Skills loaded:**
- `claude-skills-overview-2026`
- `claude-plugins-reference-2026`
- `claude-hooks-reference-2026`

**Domain knowledge provided by skills:**
- `claude-skills-overview-2026` — Complete reference for Claude Code skills system: frontmatter schema, YAML syntax rules, multiline indicator prohibitions, tool/allowed-tools field format, progressive disclosure patterns, `name:` field requirements per agentskills.io spec
- `claude-plugins-reference-2026` — Complete reference for plugin system: plugin.json schema, installation scopes, component path fields, path behavior rules, environment variables, LSP server integration, CLI commands
- `claude-hooks-reference-2026` — Complete reference for hooks system: all hook events, configuration fields, exit codes, matcher syntax

**Model:** sonnet

**Research task fitness:**

| Task | Fitness | Evidence |
|------|---------|----------|
| Search existing plugins/skills for similar functionality | HIGH — Phase 1 "Discovery" workflow explicitly Globs for all capability files | Lines 32-44: Verify plugin.json, Glob for `skills/*/SKILL.md`, `commands/*.md`, `agents/*.md` |
| Gather domain knowledge about plugin structure, conventions, frontmatter rules | HIGH — three reference skills cover all schema, YAML rules, and conventions | Loads 3 reference skills; Rules 11-12 require consulting loaded skills before making technical claims |
| Identify architecture patterns from well-structured plugins | HIGH — Phase 3 Reference File Audit + Phase 8 Cross-Reference Analysis build structural maps | Lines 100-158 for reference classification; lines 280-310 for cross-reference graph |

**Assessment:** plugin-assessor is the strongest candidate. It performs structural discovery as its primary workflow, carries all three authoritative domain-knowledge skills, and has explicit rules requiring it to cite those skills rather than guess.

---

### 2. refactor-planner

**File:** `plugins/plugin-creator/agents/refactor-planner.md`

**Skills loaded:** none (no skills field in frontmatter)

**Domain knowledge provided by skills:** None — no skills field. All domain knowledge must come from the model's training data.

**Model:** sonnet

**Research task fitness:**

| Task | Fitness | Evidence |
|------|---------|----------|
| Search existing plugins/skills for similar functionality | MEDIUM — Discovery Phase reads plugin.json and Globs component directories | Lines 21-26 perform file discovery |
| Gather domain knowledge about plugin structure, conventions, frontmatter rules | LOW — no reference skills loaded; relies on training data | No skills field; no authoritative source for schema claims |
| Identify architecture patterns from well-structured plugins | MEDIUM — Skill Analysis phase reads SKILL.md files and evaluates domains | Lines 27-35 evaluate each skill for domains and split candidates |

**Assessment:** Useful for planning after discovery is complete, but carries no authoritative schema references. Would need to guess at conventions rather than cite them.

---

### 3. subagent-refactorer

**File:** `plugins/plugin-creator/agents/subagent-refactorer.md`

**Skills loaded:**
- `write-frontmatter-description`

**Domain knowledge provided by skills:**
- `write-frontmatter-description` — Formatting rules for frontmatter description fields: single-line only, no colons, front-load critical info, 1024-char limit

**Model:** sonnet

**Tools:** Extensive — includes WebFetch, WebSearch, MCP Ref tools, context7, exa, github, sequential-thinking, episodic-memory

**Research task fitness:**

| Task | Fitness | Evidence |
|------|---------|----------|
| Search existing plugins/skills for similar functionality | LOW — designed for prompt engineering research, not codebase discovery | Phase 1 focuses on Anthropic documentation URLs, not local filesystem search |
| Gather domain knowledge about plugin structure, conventions, frontmatter rules | MEDIUM — can research Anthropic docs via MCP tools but lacks the three plugin reference skills | Phase 1 lists specific Anthropic doc URLs to consult |
| Identify architecture patterns from well-structured plugins | LOW — targets agent refactoring patterns, not plugin structural patterns | Focus is prompt engineering methodology, not plugin architecture |

**Assessment:** Wrong specialization for Phase 1 plugin research. Strong for improving agent prompts; weak for discovering local plugin conventions.

---

### 4. agent-creator

**File:** `plugins/plugin-creator/agents/agent-creator.md`

**Skills loaded:**
- `plugin-creator:claude-plugins-reference-2026`
- `plugin-creator:claude-hooks-reference-2026`
- `plugin-creator:claude-skills-overview-2026`
- `plugin-creator:agent-creator`

**Domain knowledge provided by skills:**
- Three authoritative references identical to plugin-assessor (skills, plugins, hooks)
- Plus `agent-creator` skill (creation workflow guidance)

**Model:** sonnet

**Tools:** Read, Write, Edit, Grep, Glob, Bash, Task

**Research task fitness:**

| Task | Fitness | Evidence |
|------|---------|----------|
| Search existing plugins/skills for similar functionality | HIGH — Phase 1 Discovery explicitly Globs `agents/*.md` across project and plugin directories | Lines 33-36: `Glob("agents/*.md", ".claude/")` and `Glob("plugins/*/agents/*.md")` |
| Gather domain knowledge about plugin structure, conventions, frontmatter rules | HIGH — three reference skills provide complete authoritative coverage | Frontmatter Constraints section cites exact schema requirements backed by loaded skills |
| Identify architecture patterns from well-structured plugins | HIGH — discovery reads existing agents and identifies patterns for template selection | Phase 3 Template Selection reads existing agents as templates |

**Assessment:** Strong research candidate. Carries the same three reference skills as plugin-assessor. Discovery workflow is narrower (agents only vs. all components), but it explicitly reads existing files to identify patterns.

---

### 5. refactor-validator

**File:** `plugins/plugin-creator/agents/refactor-validator.md`

**Skills loaded:** none (no skills field in frontmatter)

**Domain knowledge provided by skills:** None.

**Model:** sonnet

**Research task fitness:**

| Task | Fitness | Evidence |
|------|---------|----------|
| Search existing plugins/skills for similar functionality | LOW — validates completed refactoring, does not discover existing work | Workflow is task completion check + structure validation, not discovery |
| Gather domain knowledge about plugin structure, conventions, frontmatter rules | LOW — no reference skills; validates against hardcoded criteria | Lines 60-86 list criteria but without authoritative skill backing |
| Identify architecture patterns from well-structured plugins | LOW — checks for regressions, does not identify patterns | Phase 4 regression detection, not pattern research |

**Assessment:** Post-implementation validator. Not suited for Phase 1 research.

---

### 6. hook-creator

**File:** `plugins/plugin-creator/agents/hook-creator.md`

**Skills loaded:**
- `plugin-creator:claude-hooks-reference-2026`
- `plugin-creator:hooks-core-reference`
- `plugin-creator:hooks-io-api`
- `plugin-creator:hooks-patterns`
- `plugin-creator:hook-creator`

**Domain knowledge provided by skills:**
- `claude-hooks-reference-2026` — Complete hooks reference
- `hooks-core-reference` — Core hooks documentation
- `hooks-io-api` — Hooks I/O API spec
- `hooks-patterns` — Reusable hook patterns
- `hook-creator` — Creation workflow

**Model:** sonnet

**Research task fitness:**

| Task | Fitness | Evidence |
|------|---------|----------|
| Search existing plugins/skills for similar functionality | LOW — domain is hook scripts, not general plugin discovery | Workflow reads hooks.json and hook scripts, not skills/agents |
| Gather domain knowledge about plugin structure, conventions, frontmatter rules | LOW — authoritative only on hook domain; no plugin/skill reference skills | No claude-skills-overview or claude-plugins-reference loaded |
| Identify architecture patterns from well-structured plugins | LOW — hook engineering specialist, not plugin architecture analyst | Pattern library covers hook patterns only |

**Assessment:** Excellent for hook work; wrong specialization for Phase 1 plugin research.

---

### 7. refactor-executor

**File:** `plugins/plugin-creator/agents/refactor-executor.md`

**Skills loaded:** none (no skills field in frontmatter)

**Domain knowledge provided by skills:** None.

**Model:** sonnet

**Research task fitness:**

All three tasks: LOW — this is an execution orchestrator that runs pre-planned tasks. It does not perform discovery, gather domain knowledge, or identify patterns.

**Assessment:** Implementation agent; does not serve research purposes.

---

### 8. contextual-ai-documentation-optimizer

**File:** `plugins/plugin-creator/agents/contextual-ai-documentation-optimizer.md`

**Skills loaded:**
- `prompt-optimization-claude-45`
- `write-frontmatter-description`
- `subagent-contract`
- `plugin-creator:audit-skill-completeness`

**Domain knowledge provided by skills:**
- `prompt-optimization-claude-45` — Prompt optimization principles for Claude models
- `write-frontmatter-description` — Frontmatter description formatting rules
- `subagent-contract` — Subagent behavioral discipline
- `audit-skill-completeness` — 8-category skill quality evaluation

**Model:** sonnet

**Research task fitness:**

| Task | Fitness | Evidence |
|------|---------|----------|
| Search existing plugins/skills for similar functionality | LOW — designed for documentation optimization, not codebase discovery | Workflow is RT-ICA → analyze → optimize → verify |
| Gather domain knowledge about plugin structure, conventions, frontmatter rules | MEDIUM — knows frontmatter formatting rules and skill completeness criteria, but lacks the three full reference skills | Has write-frontmatter-description and audit-skill-completeness; lacks claude-skills-overview-2026 and claude-plugins-reference-2026 |
| Identify architecture patterns from well-structured plugins | LOW — evaluates documentation quality of a given file, does not survey the plugin landscape | File-type-specific strategies are for optimizing content, not discovering patterns |

**Assessment:** Excellent for optimizing documentation after research is complete; not suited for research itself.

---

## Comparison Matrix

| Agent | Skills Covering Schema | Discovery Workflow | Pattern Identification | Authoritative Claims Rule |
|-------|----------------------|--------------------|----------------------|--------------------------|
| plugin-assessor | 3/3 (skills, plugins, hooks) | Full multi-component Glob | Phase 8 Cross-Reference Analysis | Explicit (Rules 11-13) |
| agent-creator | 3/3 (skills, plugins, hooks) | Agent-focused Glob | Template pattern reading | Implicit via loaded skills |
| refactor-planner | 0/3 | Reads plugin.json + Globs | Skill domain analysis | None |
| subagent-refactorer | 0/3 (description rules only) | Anthropic docs only | None | Research-first mandate |
| contextual-ai-documentation-optimizer | 0/3 (formatting + completeness only) | None | None | CoVe post-check |
| hook-creator | 1/3 (hooks only) | Hook scripts only | Hook patterns only | None |
| refactor-validator | 0/3 | None | None | None |
| refactor-executor | 0/3 | None | None | None |

---

## Recommendation

**Primary recommendation: `plugin-creator:plugin-assessor`**

**Rationale:**
- Loads all three authoritative reference skills (`claude-skills-overview-2026`, `claude-plugins-reference-2026`, `claude-hooks-reference-2026`) providing complete schema, convention, and frontmatter knowledge
- Phase 1 Discovery workflow explicitly designed for comprehensive filesystem scanning across all component types (skills, commands, agents, config files)
- Phase 8 Cross-Reference Analysis builds a structural link graph — directly addresses "identify architecture patterns"
- Rules 11-13 explicitly require consulting loaded skills before making claims and citing authoritative sources — reduces hallucination risk on schema questions
- Report output documents findings with file:line citations, making research artifacts usable in subsequent phases

**Secondary option: `plugin-creator:agent-creator`**

Use agent-creator when Phase 1 research scope is limited to agents specifically (not full plugin components). It carries the same three reference skills and explicitly reads existing agents to identify templates and patterns.

**Do not use for Phase 1:**
- `refactor-planner` — no reference skills; guesses at conventions
- `subagent-refactorer` — wrong domain (prompt engineering research, not local discovery)
- `contextual-ai-documentation-optimizer` — optimization specialist, not researcher
- `hook-creator` — hook-domain only
- `refactor-validator` — post-implementation validator
- `refactor-executor` — execution orchestrator

---

## Delegation Format for SKILL.md Workflows

When referencing Phase 1 research in a plugin-creation SKILL.md workflow:

```text
1. Task is domain research and code discovery with subagent_type="plugin-creator:plugin-assessor"
   Context to include in the prompt: target plugin directory path, list of component types to survey,
     any specific conventions or patterns to look for
   Output: .claude/audits/plugin-assessment-{plugin-slug}.md — structural assessment with
     file:line references for all discovered patterns, conventions, and similar functionality
```

If agent-scope is limited to agents only:

```text
1. Task is agent pattern discovery with subagent_type="plugin-creator:agent-creator"
   Context to include in the prompt: plugin directory path, description of agent functionality needed
   Output: summary of existing agent patterns and recommended template approach
```
