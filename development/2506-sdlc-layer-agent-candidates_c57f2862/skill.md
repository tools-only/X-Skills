# SDLC Layer Separation Architecture — Agent Candidates

**Source**: Full read of all 15 agents in `.claude/agents/`
**Plan Context**: Layer 0 (SDLC-Agnostic), Layer 1 (Language-Specific), Layer 2 (Stack/Goal-Specific), ARL Meta-Layer
**Generated**: 2026-02-23

---

## Layer 0 (SDLC-Agnostic)

### backlog-item-groomer

- **Type**: agent
- **Path**: `.claude/agents/backlog-item-groomer.md`
- **Proposed Layer**: L0
- **Relevance**: Produces context manifest with RT-ICA assessment, dependency graph, and blockers. Core to SAM planning pipeline and human touchpoints before task delegation.
- **Key content**:
  - RT-ICA procedure: Goal statement → Reverse prerequisites → Availability check (AVAILABLE/DERIVABLE/MISSING) → BLOCKED/APPROVED decision
  - Discovers supporting skills, related agents, prior work, dependencies, blockers
  - Output format: Context Manifest with RT-ICA Summary, Supporting Skills, Related Agents, Prior Work, Dependencies, Blockers, Suggested First Steps
  - Uses Glob/Grep/Read only; no write tools; orchestrator owns commits

### fact-checker

- **Type**: agent
- **Path**: `.claude/agents/fact-checker.md`
- **Proposed Layer**: L0
- **Relevance**: Implements verification protocol and evidence discipline. Returns VERIFIED/REFUTED/INCONCLUSIVE with citations. Training data recall explicitly rejected.
- **Key content**:
  - Mandatory tool usage: WebFetch, WebSearch, gh, CLI — must use at least one before verdict
  - Chain of Verification (CoVe): 2–3 falsification questions, cross-check with different source/method
  - Structured verdict format: CLAIM, VERDICT, EVIDENCE, CROSS_CHECK, EXPLANATION, CITATION
  - Prohibited: VERIFIED/REFUTED without tool evidence; phrases "I know", "from my training"
  - Single-claim scope; does not update backlog or fix docs

### topic-specialist

- **Type**: agent
- **Path**: `.claude/agents/topic-specialist.md`
- **Proposed Layer**: L0
- **Relevance**: Embodies evidence discipline and CoVe for domain research. Primary sources only; can update/create skills with verified findings.
- **Key content**:
  - Invocation format: TOPIC, SKILLS, QUESTION, CONDITIONS, OUTPUT (answer only | update skill | create skill)
  - Source priority: GitHub source → README/docs → issues → official docs → WebSearch
  - Minimum research checklist: repo URL, README, main entry point, GitHub issues
  - CoVe step: 2–3 falsification questions, cross-check, revise if discrepancy
  - Evidence rules: Never state fact without URL + access date; no "I know", "typically"; quote directly
  - Boundaries: No training-data-only responses; no commits; modifies only skills/research

### context-gathering

- **Type**: agent
- **Path**: `.claude/agents/context-gathering.md`
- **Proposed Layer**: L0
- **Relevance**: Creates context manifest for new tasks. Human touchpoint for task prep; ensures developer has everything needed before implementation.
- **Key content**:
  - Process: Understand task → Research everything (spare no tokens) → Write narrative context manifest
  - Output: "How It Currently Works" (verbose narrative), "For New Feature Implementation", "Technical Reference Details"
  - CRITICAL: May ONLY use Edit on the task file; forbidden from editing other files
  - Self-verification: "Could someone implement with ONLY my context manifest?"
  - Architecture patterns, data access, code organization, business logic extraction

### context-refinement

- **Type**: agent
- **Path**: `.claude/agents/context-refinement.md`
- **Proposed Layer**: L0
- **Relevance**: Updates context manifest with discoveries from work session. Only updates if drift or new discoveries found. Accumulates institutional knowledge.
- **Key content**:
  - Process: Read transcript → Analyze for drift/discoveries → Decision (no update vs. update)
  - Drift types: Component behavior different, gotchas, hidden dependencies, wrong assumptions, environmental requirements
  - Update format: "Discovered During Implementation" with narrative + Updated Technical Details
  - Qualifies: Undocumented interactions, incorrect assumptions, hidden side effects; NOT minor typos, temporary workarounds

### logging

- **Type**: agent
- **Path**: `.claude/agents/logging.md`
- **Proposed Layer**: L0
- **Relevance**: Consolidates work logs into task Work Log section. Artifact conventions for task format; invoked during context compaction or task completion.
- **Key content**:
  - Responsibilities: Read target file + transcript → Assess cleanup → Remove outdated → Update existing → Add new → Chronological order
  - Work Log format: Date → Completed, Decisions, Discovered, Next Steps
  - Rules: Cleanup first, chronological integrity, consolidation, clarity
  - CRITICAL: May ONLY edit the specific task file; never sessions/state/, current-task.json, system state
  - Scope: Work Log, Success Criteria, Next Steps, Context Manifest

### code-review

- **Type**: agent
- **Path**: `.claude/agents/code-review.md`
- **Proposed Layer**: L0
- **Relevance**: Quality gate for code changes. Invoked explicitly or by protocol. Reviews for security, bugs, performance, project patterns. LLM-slop detection.
- **Key content**:
  - Use ONLY when explicitly requested or by protocol in sessions/protocols/
  - LLM slop focus: Reimplementing existing, failing norms, junk patterns, placeholders, hallucinated defaults, duplicate env vars
  - Categorization: Critical (blocks deployment), Warning (should address), Suggestion (consider)
  - "Keep it real": Consider actual risk; don't concern troll
  - Output: Summary, Critical/Warning/Suggestion sections, Patterns Followed, Overall Assessment

---

## Layer 1 (Language-Specific)

### javascript-pro

- **Type**: agent
- **Path**: `.claude/agents/javascript-pro.md`
- **Proposed Layer**: L1
- **Relevance**: Language-specific specialist for ES2023+, Node.js. Abstract role with quality gates and project detection.
- **Key content**:
  - Operational protocol: Read package.json, build config, module setup → Analyze patterns → Implement → Verify
  - Quality checklist: ESLint, Prettier, tests >85%, JSDoc, bundle size, error handling
  - Language standards: ES2023+ features, async patterns, module design, functional/OOP patterns
  - Project detection: package.json, webpack/rollup/vite, eslint.config, test config
  - Persistent agent memory for project-scope patterns

### typescript-pro

- **Type**: agent
- **Path**: `.claude/agents/typescript-pro.md`
- **Proposed Layer**: L1
- **Relevance**: Language-specific specialist for TypeScript 5.0+. Abstract role with quality gates and project detection.
- **Key content**:
  - Initialization protocol: Read tsconfig.json, package.json, build configs → Assess type patterns → Identify framework
  - TypeScript checklist: Strict mode, no any, 100% type coverage, type-only imports, source maps
  - Advanced patterns: Conditional types, mapped types, branded types, Result types, discriminated unions
  - Quality verification: tsc --noEmit, lint, type coverage, tests, bundle analysis
  - Project detection: tsconfig.json, framework (React/Vue/Node), lint config

### c-systems-programmer

- **Type**: agent
- **Path**: `.claude/agents/c-systems-programmer.md`
- **Proposed Layer**: L1
- **Relevance**: Language-specific specialist for C (C99/C11), systems programming, embedded. Abstract role with quality gates.
- **Key content**:
  - Technical stack: C99/C11, POSIX, GCC/Clang, Valgrind, GDB, gprof, perf
  - Review checklist: Memory management, pointer safety, error handling, thread safety, performance, code quality
  - Quality gates: -Wall -Wextra, Valgrind --leak-check=full, clang-tidy
  - Design patterns: Memory pool, error handling in syscalls, thread-safe singleton, lock-free structures

---

## Layer 2 (Stack/Goal-Specific)

### plugin-assessor

- **Type**: agent
- **Path**: `.claude/agents/plugin-assessor.md`
- **Proposed Layer**: L2
- **Relevance**: Stack-specific (Claude Code plugins). Validates schema, frontmatter, reference organization. Architecture patterns for plugin structure.
- **Key content**:
  - Phases: Discovery → Manifest Validation → Skills Analysis → Commands → Agents → Hooks → MCP → Cross-Reference → Enhancement
  - Schema validation: plugin.json, skill frontmatter, command/agent frontmatter
  - Orphan classification: Linked, Orphaned-New Content, Duplicate, Notes, Examples
  - Scoring: Structural validity, manifest completeness, frontmatter, description quality, reference organization
  - Skills: claude-skills-overview-2026, claude-plugins-reference-2026, claude-hooks-reference-2026

### plugin-docs-writer

- **Type**: agent
- **Path**: `.claude/agents/plugin-docs-writer.md`
- **Proposed Layer**: L2
- **Relevance**: Stack-specific. Generates user-facing README for Claude Code plugins. Translates AI-facing content to human outcomes.
- **Key content**:
  - Research process: plugin.json → each skill SKILL.md + references → commands → hooks → MCP
  - Translation: Skill behavior → "Claude will..."; Commands → syntax + args; Hooks → "Automatically..."
  - Banned terms: ROLE_TYPE, orchestrator, frontmatter, scientific method, permissionMode
  - README structure: Why Install, What You Get, Commands, Claude Improvements, Automatic Behaviors, Installation, Usage, Example

### research-curator

- **Type**: agent
- **Path**: `.claude/agents/research-curator.md`
- **Proposed Layer**: L2
- **Relevance**: Creates structured research entries for tools/libraries. Stack research into knowledge base. Operates standalone or orchestrated.
- **Key content**:
  - Modes: Default (new URL), --rerun (re-research), --fix (validation issues)
  - Tools: MCP Ref, exa, gh API for repo metadata
  - Categories: research-agent-patterns, mcp-ecosystem, agent-frameworks, etc.
  - Entry template: 10 sections, no placeholders, freshness tracking (3 months)
  - Boundaries: Does not update README, commit, or modify outside research/

### research-context-agent

- **Type**: agent
- **Path**: `.claude/agents/research-context-agent.md`
- **Proposed Layer**: L2
- **Relevance**: Cross-references research with skills, agents, hooks, commands. Discovers integration opportunities. Connective tissue between research and capabilities.
- **Key content**:
  - Three phases: Absorb (extract from research) → Search & Match (5 dimensions) → Append (Integration Opportunities section)
  - Dimensions: Enhance skills, enhance agents, enhance hooks, enhance commands, new skill candidate, new MCP candidate
  - Output: Enhances Existing table, New Skill Candidates, New MCP Server Candidates, Cross-References
  - Rules: Concrete over vague, skip empty sections, no false positives, preserve content, idempotent

---

## ARL Meta-Layer

### doc-drift-auditor

- **Type**: agent
- **Path**: `.claude/agents/doc-drift-auditor.md`
- **Proposed Layer**: ARL (Identify)
- **Relevance**: Identifies deviations between documentation and implementation. Git forensics + code analysis. Evidence-based drift reporting.
- **Key content**:
  - Process: Repository discovery → Git timeline → Implementation analysis → Documentation claims → Drift detection → Report
  - Drift categories: Implemented but undocumented, Documented but unimplemented, Outdated, Mismatched details
  - Evidence: file:line, commit SHA, quoted claims, code reality
  - Output: DOCUMENTATION_DRIFT_AUDIT.md with executive summary, timeline, findings by category
  - Boundaries: Audit only; no modifications; no assumptions; no training-data reliance

### code-review (ARL overlap)

- **Type**: agent
- **Path**: `.claude/agents/code-review.md`
- **Proposed Layer**: ARL (Identify)
- **Relevance**: LLM-slop and hallucination detection. Identifies reimplemented scaffolding, junk patterns, hallucinated defaults.
- **Key content**:
  - Primary focus: "Some or all of the code you are reviewing was generated by an LLM"
  - Hallucination patterns: Placeholders, TODOs, imagined defaults, duplicate env vars
  - Complements doc-drift-auditor for code vs. spec drift

### context-refinement (ARL overlap)

- **Type**: agent
- **Path**: `.claude/agents/context-refinement.md`
- **Proposed Layer**: ARL (Accumulate)
- **Relevance**: Accumulates local domain knowledge discovered during implementation. Reduces staleness of context manifest.
- **Key content**:
  - Captures: Undocumented interactions, incorrect assumptions, hidden side effects
  - "Guardian of institutional knowledge" — helps future developers avoid same surprises

### logging (ARL overlap)

- **Type**: agent
- **Path**: `.claude/agents/logging.md`
- **Proposed Layer**: ARL (Observe)
- **Relevance**: Observes and consolidates work session output. Maintains clean task state for future sessions.
- **Key content**:
  - Consolidates Completed, Decisions, Discovered, Next Steps
  - Ensures task file reflects present state for kaizen/session-historian consumption

### topic-specialist (ARL overlap)

- **Type**: agent
- **Path**: `.claude/agents/topic-specialist.md`
- **Proposed Layer**: ARL (Improve)
- **Relevance**: Can update or create skills with verified knowledge. Feeds SAM + connective tissues with primary-source-verified content.
- **Key content**:
  - OUTPUT: "answer + update skill" or "answer + create skill"
  - Update: Append to SKILL.md with citations; do not remove existing
  - Create: Invoke skill-creator, add-doc-updater, populate with verified findings

### research-context-agent (ARL overlap)

- **Type**: agent
- **Path**: `.claude/agents/research-context-agent.md`
- **Proposed Layer**: ARL (Accumulate)
- **Relevance**: Accumulates integration opportunities from research into structured Cross-References. Connective tissue between research KB and skills/agents.
- **Key content**:
  - Appends "Integration Opportunities" to research files
  - Cross-references between research files (e.g., logfire + tensorzero)

---

## Summary by Layer

| Layer | Agents |
|-------|--------|
| **L0** | backlog-item-groomer, fact-checker, topic-specialist, context-gathering, context-refinement, logging, code-review |
| **L1** | javascript-pro, typescript-pro, c-systems-programmer |
| **L2** | plugin-assessor, plugin-docs-writer, research-curator, research-context-agent |
| **ARL** | doc-drift-auditor (Identify), code-review (Identify), context-refinement (Accumulate), logging (Observe), topic-specialist (Improve), research-context-agent (Accumulate) |

---

## Cross-Cutting Patterns

1. **Evidence discipline**: fact-checker, topic-specialist — primary sources only, CoVe, citations
2. **RT-ICA**: backlog-item-groomer — pre-planning gate, AVAILABLE/DERIVABLE/MISSING
3. **Context manifest**: backlog-item-groomer, context-gathering, context-refinement — task format artifact
4. **Orchestrator boundaries**: Most agents defer commits, README updates, batch coordination to orchestrator
5. **Project detection**: L1 agents (javascript-pro, typescript-pro, c-systems-programmer) read config files before acting
6. **Stack-specific skills**: plugin-assessor, plugin-docs-writer load claude-skills-overview-2026, claude-plugins-reference-2026, claude-hooks-reference-2026
