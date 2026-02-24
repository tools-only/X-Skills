# SDLC Layer 2 Integration Suggestions

**Generated**: 2026-02-23  
**Plan Context**: SDLC Layer Separation Architecture — Layer 2 = stack profiles, architecture patterns, toolchain config, reference examples, stack research  
**Source**: Full read of each Layer 2 item from sdlc-layer-candidates-master.md

---

## How to Use This Document

For each Layer 2 item:
- **Integration**: Where in the plan this item should be placed
- **Amendment**: Proposed change to the plan
- **Nugget**: Detail for other steps (e.g., implementation, schema design)

---

## 1. the-rewrite-room (Plugin)

**Path**: `plugins/the-rewrite-room/skills/the-rewrite-room/SKILL.md`

### Extracted Content

Workflow router for documentation and authoring. Maps task intent to 7 canonical workflows: drift-audit, documentation-sync, authoring, prompt-optimization, summarization, formatting-validation, research-utilities. Uses STATUS block output contract (DONE/BLOCKED/FAILED). Adapter shims normalize legacy component output. Router script classifies tasks; link_checker, file_metrics support validation.

### Integration

**Where**: Layer 2 schema — "Documentation & Authoring Stack Profile"

Add a stack profile entry for documentation workflows that routes through the-rewrite-room. The router taxonomy (drift-audit → doc-drift-auditor, prompt-optimization → optimize-claude-md, etc.) is a reference architecture pattern for intent-based workflow routing.

### Amendment

Add to Layer 2 schema:

```yaml
stack_profiles:
  documentation-authoring:
    router: the-rewrite-room
    workflows: [drift-audit, documentation-sync, authoring, prompt-optimization, summarization, formatting-validation, research-utilities]
    output_contract: STATUS block (DONE|BLOCKED|FAILED)
    canonical_agents: [doc-drift-auditor, contextual-ai-documentation-optimizer]
    canonical_skills: [optimize-claude-md, summarizer, add-doc-updater, plugin-creator:lint]
```

### Nugget

The STATUS block format (STATUS, SUMMARY, ARTIFACTS, VALIDATION, NOTES) is a reusable output contract for Layer 2 workflows. Consider standardizing this across stack profiles.

---

## 2. the-rewrite-room/workflows (Plugin)

**Paths**: `plugins/the-rewrite-room/skills/the-rewrite-room/workflows/*.md` (drift-audit, documentation-sync, authoring, prompt-optimization, summarization, formatting-validation, research-utilities)

### Extracted Content

Each workflow has YAML frontmatter: workflow name, canonical_agent/skill, canonical_path, version, output_contract. Entrypoint contracts define required/optional inputs. Steps, validation gates (HARD STOP, SOFT STOP), and output contracts are specified. Adapters map legacy output to STATUS block.

### Integration

**Where**: Layer 2 reference examples — "Workflow Definition Schema"

Use these as the canonical format for Layer 2 workflow definitions. The frontmatter schema (workflow, canonical_*, output_contract) should inform the Layer 2 workflow registry schema.

### Amendment

Define a Layer 2 workflow schema derived from the-rewrite-room frontmatter:

```yaml
workflow_schema:
  required: [workflow, version, output_contract]
  optional: [canonical_agent, canonical_skill, canonical_path, knowledge_reference]
  validation_gates: [HARD_STOP, SOFT_STOP]
```

### Nugget

The formatting-validation workflow references both plugin-creator:lint and holistic-linting as complementary validators. Layer 2 toolchain config should support multi-validator composition.

---

## 3. summarizer (Plugin)

**Path**: `plugins/summarizer/skills/summarizer/SKILL.md`

### Extracted Content

Format routing (structured, bullets, tldr, json, table, outline). Decision tree: file → file-summarization; URL → url-summarization; image → image-summarization; multi-source → team coordination (3+ sources). Seven fidelity rules (read before summarizing, extract before abstracting, preserve counts, etc.). Structured output with YAML frontmatter. Team coordination via Teammate tool for 3+ sources.

### Integration

**Where**: Layer 2 stack profile — "Summarization Stack" (subset of documentation-authoring)

The summarizer is a canonical skill for the summarization workflow. Its fidelity rules and format templates are stack-specific constraints.

### Amendment

Document summarizer as a Layer 2 component with:
- Format templates as stack config
- Fidelity rules as mandatory constraints
- Team coordination pattern (3+ sources → Teammate) as architecture pattern

### Nugget

Fidelity rules (read before summarizing, extract before abstracting, preserve counts, distinguish absence from nonexistence, no lossy re-summarization) are anti-hallucination constraints. Layer 2 schema could define a `fidelity_rules` array for summarization-type workflows.

---

## 4. holistic-linting (Plugin)

**Path**: `plugins/holistic-linting/skills/holistic-linting/SKILL.md`

### Extracted Content

Lint orchestration for orchestrators vs sub-agents. Orchestrators delegate to linting-root-cause-resolver; sub-agents run format → lint → resolve directly. Linter detection via config scan (pyproject.toml, .pre-commit-config.yaml, package.json, .clang-format, etc.). Rules knowledge base: Ruff (933 rules), MyPy, Bandit. Pre-existing issues protocol: blocking vs non-blocking, record to BACKLOG.md. Never suppress; never delete to fix.

### Integration

**Where**: Layer 2 toolchain config — "Linting Stack Profile"

Holistic-linting defines a lint orchestration pattern: detection → format → lint → resolve. It integrates with clang-format (detected via .clang-format), ruff, mypy, bandit, prek/pre-commit.

### Amendment

Add lint orchestration to Layer 2 toolchain schema:

```yaml
linting_stack:
  detection: config_scan  # pyproject.toml, .pre-commit-config.yaml, etc.
  orchestrator_behavior: delegate_to_resolver_agent
  subagent_behavior: format_then_lint_then_resolve
  rules_kb: [ruff, mypy, bandit]
  pre_existing_protocol: record_to_backlog
```

### Nugget

The holistic-linting ↔ clang-format relationship: holistic-linting detects .clang-format and includes C/C++ in its scope. Layer 2 should model toolchain composition (e.g., holistic-linting + clang-format for C++ projects).

---

## 5. gitlab-skill (Plugin)

**Path**: `plugins/gitlab-skill/skills/gitlab-skill/SKILL.md`

### Extracted Content

Four domains: CI/CD pipeline config (.gitlab-ci.yml), GLFM documentation, local pipeline testing (gitlab-ci-local), GitLab CLI (glab). validate_glfm.py for GLFM rendering. sync-gitlab-docs.py for doc updates. CRITICAL_SYNTAX_RULES for GLFM (alert blocks lowercase, collapsible single-line). glab ci lint, glab ci status (non-interactive). Documentation index from upstream sync.

### Integration

**Where**: Layer 2 stack profile — "GitLab CI/CD & GLFM Stack"

GitLab-specific stack: CI config, GLFM docs, local testing, glab CLI. Complements GitHub (gh skill) as alternative VCS/CI stack.

### Amendment

Add GitLab stack profile:

```yaml
stack_profiles:
  gitlab:
    domains: [ci_pipeline, glfm_docs, local_testing, glab_cli]
    tools: [validate_glfm.py, sync-gitlab-docs.py, gitlab-ci-local, glab]
    constraints: CRITICAL_SYNTAX_RULES
```

### Nugget

formatting-validation workflow chains to glfm-validator when GLFM content is present. Layer 2 workflow chaining should support conditional validators (e.g., if target contains GLFM → run validate_glfm.py).

---

## 6. fastmcp-creator (Plugin)

**Path**: `plugins/fastmcp-creator/skills/fastmcp-creator/SKILL.md`

### Extracted Content

MCP server creation: FastMCP 3.x (Python), TypeScript SDK. Four-phase workflow: Deep Research → Implementation → Review → Evaluations. Agent-centric design (workflow tools, not API wrappers; optimize for context; actionable errors). Provider/transform architecture, component versioning, session state, authorization. Evaluation harness for testing server quality. Standalone skill, no external dependencies.

### Integration

**Where**: Layer 2 stack profile — "MCP Server Development Stack"

FastMCP-creator is the reference stack for building MCP servers. It defines toolchain (FastMCP 3.x, TypeScript SDK), architecture patterns (provider/transform, agent-centric design), and quality gates (evaluation creation).

### Amendment

Add MCP development stack:

```yaml
stack_profiles:
  mcp_server_development:
    frameworks: [FastMCP 3.x, TypeScript SDK]
    workflow_phases: [research, implementation, review, evaluations]
    design_principles: agent_centric
    quality_gate: evaluation_harness
```

### Nugget

FastMCP-creator activates python3-development for project setup. Layer 2 stack profiles should declare dependencies on Layer 1 (language) profiles.

---

## 7. dasel (Plugin)

**Path**: `plugins/dasel/skills/dasel-reference/SKILL.md`

### Extracted Content

dasel v3 CLI for structured data (JSON, YAML, TOML, XML, CSV, HCL, INI). No -f/--file; input from stdin. Selectors, functions, modification syntax. Format conversion. Domain skills: data-exploration, data-transformation, enterprise-* (Hibernate, Spring, Tomcat, Maven). data-explorer agent.

### Integration

**Where**: Layer 2 toolchain config — "Structured Data Stack"

dasel is a stack-specific tool for query/transform of config and data files. Used by data-explorer agent, format-recipes, transformation-patterns.

### Amendment

Add structured data stack:

```yaml
toolchain_config:
  structured_data:
    primary_tool: dasel
    formats: [json, yaml, toml, xml, csv, hcl, ini]
    agents: [data-explorer, data-analyst]
    domain_skills: [data-exploration, data-transformation]
```

### Nugget

dasel v3 breaking changes (no -f flag, stdin-only) are documented. Layer 2 toolchain config should support version-specific behavior and migration notes.

---

## 8. clang-format (Plugin)

**Path**: `plugins/clang-format/skills/clang-format/SKILL.md`

### Extracted Content

C/C++ formatting. Templates (google-cpp-modified, linux-kernel, microsoft-visual-studio, modern-cpp17-20, etc.). Integration scripts: pre-commit, vimrc, emacs. Reference docs by category (alignment, breaking, braces, indentation, spacing, includes, languages, comments, advanced). Code style analysis workflow with impact scoring. Minimal-disruption requests.

### Integration

**Where**: Layer 2 toolchain config — "C/C++ Formatting Stack"

clang-format is a language-specific (C/C++) formatting tool. Composes with holistic-linting (detected via .clang-format).

### Amendment

Add C/C++ formatting to toolchain:

```yaml
toolchain_config:
  c_cpp_formatting:
    tool: clang-format
    templates: assets/configs/
    integrations: [pre-commit, vim, emacs]
    detection_marker: .clang-format
```

### Nugget

clang-format's "Analyzing Existing Code Style" workflow uses impact scoring (line changes × 10 + whitespace × 1). Layer 2 could define impact metrics for toolchain introduction.

---

## 9. xdg-base-directory (Plugin)

**Path**: `plugins/xdg-base-directory/skills/xdg-base-directory/SKILL.md`

### Extracted Content

XDG Base Directory Specification v0.8. Environment variables: XDG_CONFIG_HOME, XDG_DATA_HOME, XDG_STATE_HOME, XDG_CACHE_HOME, XDG_RUNTIME_DIR. Stdlib-only and platformdirs implementations. Directory selection guide (config vs data vs cache vs state). Anti-patterns (legacy ~/.appname, ignoring env vars, relative paths).

### Integration

**Where**: Layer 2 architecture pattern — "Cross-Platform Config Storage"

xdg-base-directory is a cross-cutting pattern for where applications store config, data, cache, state. Relevant for CLI tools, MCP servers, plugins that persist state.

### Amendment

Add architecture pattern:

```yaml
architecture_patterns:
  config_storage:
    spec: XDG Base Directory v0.8
    env_vars: [XDG_CONFIG_HOME, XDG_DATA_HOME, XDG_STATE_HOME, XDG_CACHE_HOME]
    implementations: [stdlib_only, platformdirs]
```

### Nugget

MCP servers and plugins that write config/cache should follow XDG. Layer 2 stack profiles for CLI/MCP tools should reference this pattern.

---

## 10. plugin-creator (Plugin)

**Path**: `plugins/plugin-creator/skills/plugin-creator/SKILL.md`

### Extracted Content

Orchestrates plugin creation via delegation. Phases: RT-ICA → Discussion → Research (4-way parallel) → Design → Validation → Implementation. Artifact system: PROJECT.md, REQUIREMENTS.md, STATE.md, design-PLAN.md, etc. Parallel agent spawning. plugin_validator.py, create_plugin.py, fix_tool_formats.py. Refactor workflows (refactor-plugin, refactor-skill). plugin-assessor, contextual-ai-documentation-optimizer agents.

### Integration

**Where**: Layer 2 stack profile — "Claude Code Plugin Development Stack"

plugin-creator is the meta-stack for building plugins. It defines the plugin lifecycle, validation toolchain, and refactoring workflows.

### Amendment

Add plugin development stack:

```yaml
stack_profiles:
  plugin_development:
    orchestrator: plugin-creator
    phases: [rt_ica, discussion, research, design, validation, implementation]
    tools: [plugin_validator.py, create_plugin.py, fix_tool_formats.py]
    refactor_workflows: [refactor-plugin, refactor-skill]
    agents: [plugin-assessor, refactor-planner, refactor-executor]
```

### Nugget

plugin-creator's artifact system (STATE.md, design-PLAN.md) aligns with SAM methodology. Layer 2 should reference SAM artifact conventions for plan/state persistence.

---

## 11. development-harness agents (Plugin)

**Paths**: `plugins/development-harness/agents/*.md` (swarm-task-planner, plan-validator, feature-researcher, codebase-analyzer, doc-drift-auditor, etc.)

### Extracted Content

swarm-task-planner: CLEAR+CoVe task design, dependency-based plans, parallelization markers, sync checkpoints. plan-validator: 9 validation dimensions (requirement coverage, task completeness, dependency correctness, role-agent match, I/O validity, artifact wiring, testability, scope sanity, boundary compliance). feature-researcher: WHO/WHAT/WHEN/WHY discovery, feature-context-{slug}.md. codebase-analyzer: PATTERNS.md, ARCHITECTURE.md, TESTING.md, CONVENTIONS.md, CONCERNS.md.

### Integration

**Where**: Layer 2 reference examples — "Development Harness Agent Archetypes"

These agents implement the SAM 7-stage pipeline. They are stack-specific in that they resolve roles via language manifest (Layer 1). swarm-task-planner, plan-validator, feature-researcher, codebase-analyzer are planning/research archetypes.

### Amendment

Document development-harness agents as Layer 2 reference implementations:

```yaml
reference_agents:
  planning:
    swarm_task_planner: CLEAR+CoVe, dependency graphs, parallelization
    plan_validator: 9 dimensions, READY/BLOCKED
  research:
    feature_researcher: discovery docs, WHO/WHAT/WHEN/WHY
    codebase_analyzer: PATTERNS, ARCHITECTURE, TESTING, CONVENTIONS, CONCERNS
```

### Nugget

plan-validator's 9 dimensions (requirement coverage, task completeness, dependency correctness, etc.) could inform a Layer 2 "plan validation schema" for any stack's task plans.

---

## 12. gh (Skill)

**Path**: `.claude/skills/gh/SKILL.md`

### Extracted Content

GitHub CLI setup and usage. -R flag required for proxy environments. Installation script with SHA256 verification. Common commands: pr, issue, label, run, release, project. Milestones, Projects V2, label taxonomy. references: milestones.md, labels.md, projects-v2.md, issue-stories.md.

### Integration

**Where**: Layer 2 stack profile — "GitHub Project Management Stack"

gh skill is the GitHub-specific toolchain for issues, PRs, milestones, Projects V2, labels. Complements gitlab-skill for GitLab.

### Amendment

Add GitHub stack:

```yaml
stack_profiles:
  github:
    tool: gh
    constraint: -R flag in proxy environments
    domains: [issues, prs, milestones, projects_v2, labels, releases]
    references: [milestones, labels, projects-v2, issue-stories]
```

### Nugget

gh references (milestones.md, projects-v2.md) are stack research. Layer 2 schema should support a `references` array linking to KB/reference docs.

---

## 13. gh/references (Skill)

**Paths**: `.claude/skills/gh/references/milestones.md`, `labels.md`, `projects-v2.md`, `issue-stories.md`

### Extracted Content

Milestone CRUD via REST API. Label taxonomy (priority, type, status). Projects V2 GraphQL, updateProjectV2ItemFieldValue. Issue lifecycle and templates.

### Integration

**Where**: Layer 2 stack research — "GitHub API Reference"

These are reference docs that support the gh skill. Layer 2 stack research should inventory such references.

### Amendment

Layer 2 stack research schema:

```yaml
stack_research:
  github:
    references: [milestones, labels, projects-v2, issue-stories]
    api_type: [REST, GraphQL]
```

### Nugget

projects-v2.md documents updateProjectV2ItemFieldValue for status updates. Backlog item "github_project_setup.py: add Projects V2 status" depends on this. Layer 2 should track implementation gaps against reference docs.

---

## 14. create-backlog-item (Skill)

**Path**: `.claude/skills/create-backlog-item/SKILL.md`

### Extracted Content

Three modes: guided (AskUserQuestion), quick {title}, --auto {title}. Writes to BACKLOG.md, updates frontmatter counts. Duplicate detection. GitHub Issue creation for P0/P1. --auto derives from research files, skips prompts.

### Integration

**Where**: Layer 2 stack profile — "Backlog Management Stack"

create-backlog-item is part of the backlog ↔ GitHub workflow. Composes with gh, work-backlog-item, create-milestone, group-items-to-milestone.

### Amendment

Add backlog stack:

```yaml
stack_profiles:
  backlog_management:
    skills: [create-backlog-item, work-backlog-item]
    workflow: [create → groom → work → close]
    github_integration: [create_issue, milestone_assignment]
```

### Nugget

--auto mode derives from research/ — Layer 2 stack profiles could declare `research_sources` for autonomous derivation.

---

## 15. create-milestone, group-items-to-milestone, start-milestone, complete-milestone (Skills)

**Paths**: `.claude/skills/create-milestone/`, `group-items-to-milestone/`, `start-milestone/`, `complete-milestone/`

### Extracted Content

create-milestone: Creates GitHub milestone, returns number. group-items-to-milestone: Assigns BACKLOG items to milestone. start-milestone: Transitions status labels to in-progress. complete-milestone: Closes milestone, offers carry-forward.

### Integration

**Where**: Layer 2 stack profile — "Milestone Workflow Stack"

These four skills form a milestone lifecycle: create → group → start → complete.

### Amendment

```yaml
stack_profiles:
  milestone_workflow:
    sequence: [create-milestone, group-items-to-milestone, start-milestone, complete-milestone]
    artifact: milestone number (from create)
    handoff: milestone number → group-items → start → complete
```

### Nugget

complete-milestone skill has a complete-milestone SKILL.md with audit logic. Layer 2 workflow definitions should reference these as canonical implementations.

---

## 16. create-merge-request-changelog (Skill)

**Path**: `.claude/skills/create-merge-request-changelog/SKILL.md`

### Extracted Content

Analyzes git branches/diffs. Domain-based change categorization (bug fixes, enhancements, technical debt, documentation, testing). analyze_git_changes.py, fetch_gitlab_mr.py. AI analysis prompts. format_mr_description.py. Works with any git repo.

### Integration

**Where**: Layer 2 stack profile — "Changelog & Release Notes Stack"

create-merge-request-changelog is the pipeline for MR/PR descriptions and release notes. Uses analyze → AI-categorize → format.

### Amendment

```yaml
stack_profiles:
  changelog_generation:
    pipeline: [analyze_git_changes, ai_analysis, format_mr_description]
    output: domain_categorized_changelog
    input_sources: [git_branch, gitlab_mr]
```

### Nugget

daily-releases uses the same pipeline. Layer 2 should document shared pipelines (analyze → categorize → format) as reusable patterns.

---

## 17. daily-releases (Skill)

**Path**: `.claude/skills/daily-releases/SKILL.md`

### Extracted Content

Creates GitHub Releases per calendar day with commits. Uses create-merge-request-changelog pipeline (analyze_git_changes, analysis_prompts, format_mr_description). Idempotent: skips up-to-date days. list_daily_ranges.py. --start-date, --end-date, --branch, --dry-run.

### Integration

**Where**: Layer 2 stack profile — "Daily Release Automation" (extends changelog stack)

daily-releases is a specialized application of the changelog pipeline for per-day releases.

### Amendment

```yaml
stack_profiles:
  daily_releases:
    extends: changelog_generation
    trigger: calendar_day_with_commits
    idempotent: true
    scripts: [list_daily_ranges, analyze_git_changes, format_mr_description]
```

### Nugget

Layer 2 stack profiles should support `extends` for composition (daily-releases extends changelog_generation).

---

## 18. agent-browser (Skill)

**Path**: `.claude/skills/agent-browser/SKILL.md`

### Extracted Content

Browser automation CLI. Commands: open, snapshot, click, fill, type, select, wait, screenshot. Element refs (@e1, @e2) from snapshot -i. Command chaining with &&. Lock/unlock workflow for MCP. Waiting strategy: incremental waits with snapshot checks.

### Integration

**Where**: Layer 2 stack profile — "Browser Automation Stack"

agent-browser is the CLI for web interaction. Used for testing, scraping, form filling.

### Amendment

```yaml
stack_profiles:
  browser_automation:
    tool: agent-browser
    workflow: navigate → snapshot → interact → re-snapshot
    mcp_integration: cursor-ide-browser
```

### Nugget

agent-browser has MCP equivalent (cursor-ide-browser). Layer 2 should document CLI vs MCP tool alternatives for the same capability.

---

## 19. orchestrating-swarms / swarm-primitives / swarm-spawning / swarm-patterns / swarm-operations (Skills)

**Paths**: `.claude/skills/orchestrating-swarms/`, `swarm-primitives/`, `swarm-spawning/`, `swarm-patterns/`, `swarm-operations/`

### Extracted Content

Facade skill (orchestrating-swarms) loads specialist skills. Primitives: teams, teammates, tasks, inboxes. Spawning: subagents vs teammates, agent types, backends. Patterns: parallel reviews, pipeline workflows, research-then-implement. Operations: TeamCreate, SendMessage, TeamDelete, Task parameters.

### Integration

**Where**: Layer 2 architecture pattern — "Multi-Agent Swarm Orchestration"

Swarm skills define the Claude Code v2.1.45+ multi-agent architecture. Stack-specific in that they target Claude Code's swarm system.

### Amendment

```yaml
architecture_patterns:
  swarm_orchestration:
    primitives: [teams, teammates, tasks, inboxes]
    skills: [swarm-primitives, swarm-spawning, swarm-patterns, swarm-operations]
    facade: orchestrating-swarms
    stack: Claude Code v2.1.45+
```

### Nugget

summarizer uses Teammate tool for 3+ sources. Layer 2 should document which stack profiles use swarm primitives (summarizer, research-curator waves, etc.).

---

## 20. external-pattern-integrator (Skill)

**Path**: `.claude/skills/external-pattern-integrator/SKILL.md`

### Extracted Content

Integrates patterns from URLs/files into local skills/agents/plugins. Three phases: Parallel Candidate Mapping → Contextual Enhancement → Validation. Tracks source → candidates. Matches by purpose, workflow stage. Enhances without breaking workflow coherence. GSD, BMAD-METHOD integration.

### Integration

**Where**: Layer 2 architecture pattern — "External Pattern Integration"

external-pattern-integrator is a meta-pattern for evolving Layer 2 content from external sources (GSD, BMAD-METHOD, etc.).

### Amendment

```yaml
architecture_patterns:
  external_integration:
    workflow: [candidate_mapping, contextual_enhancement, validation]
    sources: [GSD, BMAD-METHOD, URLs]
    output: enhanced local skills/agents/plugins
```

### Nugget

Layer 2 plan should include a step: "Run external-pattern-integrator on GSD/BMAD-METHOD to discover integration opportunities."

---

## 21. universe-of-thoughts (Skill)

**Path**: `.claude/skills/universe-of-thoughts/SKILL.md`

### Extracted Content

Creative reasoning for ill-defined problems. Three paradigms: Combinational (novel combinations), Exploratory (expand boundaries), Transformative (alter constraints). When: ambiguous goals, vast solution space, innovation needed. When NOT: well-defined, mathematical, convergent. Based on Boden's creativity theory.

### Integration

**Where**: Layer 2 architecture pattern — "Creative Reasoning Stack"

universe-of-thoughts is a reasoning framework for ill-defined problems. Complements CoVe (factual) and scientific-thinking (hypothesis-driven).

### Amendment

```yaml
architecture_patterns:
  creative_reasoning:
    framework: universe-of-thoughts
    paradigms: [combinational, exploratory, transformative]
    triggers: [ambiguous_goals, vast_solution_space, innovation_required]
    excludes: [well_defined, mathematical, convergent]
```

### Nugget

Layer 2 plan design (stack profiles, schema) may benefit from UoT when goals are ambiguous. Document when to use UoT vs RT-ICA vs CoVe in plan design.

---

## 22. plugin-assessor (Agent)

**Path**: `.claude/agents/plugin-assessor.md`

### Extracted Content

Analyzes plugins for structure, frontmatter, schema compliance, quality. Comprehensive reference file audit, orphan detection. Cross-reference validation. Generates assessment reports. Loads claude-skills-overview-2026, claude-plugins-reference-2026, claude-hooks-reference-2026.

### Integration

**Where**: Layer 2 reference agent — "Plugin Quality Assessment"

plugin-assessor validates plugin structure. Used by plugin-creator refactor workflow and pre-marketplace review.

### Amendment

```yaml
reference_agents:
  plugin_quality:
    plugin_assessor: [structure, frontmatter, orphans, cross_refs]
    consumed_by: [refactor-plugin, marketplace_prep]
```

### Nugget

plugin-assessor's orphan classification could inform Layer 2 schema validation (detect orphaned stack profile entries).

---

## 23. plugin-docs-writer (Agent)

**Path**: `.claude/agents/plugin-docs-writer.md`

### Extracted Content

Generates user-facing README.md for plugins. Researches capabilities, translates AI-facing content to human terms. Documents skills (behavioral changes → outcomes), commands (syntax, args), hooks (automation). Translation examples.

### Integration

**Where**: Layer 2 reference agent — "Documentation Generation"

plugin-docs-writer produces user-facing docs from AI-facing specs. Relevant when Layer 2 plan produces user-facing documentation.

### Amendment

```yaml
reference_agents:
  documentation:
    plugin_docs_writer: [README generation, skill→outcome translation]
```

### Nugget

Layer 2 stack profiles are AI-facing. If producing user docs (e.g., "SDLC Layer Guide for Contributors"), plugin-docs-writer pattern applies: translate stack profile → human outcomes.

---

## 24. research-curator (Agent)

**Path**: `.claude/agents/research-curator.md`

### Extracted Content

Creates research entries in ./research/. Modes: single URL, --batch, --rerun, --validate. Gathers via gh API, MCP Ref, mcp exa. Category selection. Entry template. Freshness tracking. Works standalone or orchestrated by /research-curator skill.

### Integration

**Where**: Layer 2 stack research — "Research KB Population"

research-curator populates research/ which feeds Layer 2 stack research. create-backlog-item --auto derives from research/.

### Amendment

```yaml
stack_research:
  population: research-curator agent
  location: ./research/
  format: entry_template
  freshness: tracking_section
```

### Nugget

Layer 2 stack research entries could follow research/ entry template. research-curator could create entries for new stack profiles (e.g., "research/sdlc-stacks/gitlab-ci.md").

---

## 25. research-context-agent (Agent)

**Path**: `.claude/agents/research-context-agent.md`

### Extracted Content

Cross-references research files with skills, agents, hooks, plugins. Three phases: Absorb → Search & Match → Append. Five dimensions: enhance skills, enhance agents, enhance hooks, enhance commands, new skill/MCP candidates. Appends "Integration Opportunities" section to research files.

### Integration

**Where**: Layer 2 architecture pattern — "Research-to-Stack Integration"

research-context-agent discovers integration opportunities. Run after adding Layer 2 stack profiles to find enhancement targets.

### Amendment

```yaml
architecture_patterns:
  research_integration:
    agent: research-context-agent
    output: Integration Opportunities section
    dimensions: [skills, agents, hooks, commands, new_skill, new_mcp]
```

### Nugget

After defining Layer 2 schema, run research-context-agent on plan doc or new stack research entries to discover integration opportunities with existing skills/agents.

---

## Summary: Plan Amendments

| Amendment | Description |
|-----------|-------------|
| **Stack profiles** | Add: documentation-authoring, gitlab, mcp_server_development, structured_data, c_cpp_formatting, plugin_development, github, backlog_management, milestone_workflow, changelog_generation, daily_releases, browser_automation |
| **Architecture patterns** | Add: config_storage (XDG), swarm_orchestration, external_integration, creative_reasoning, research_integration |
| **Toolchain config** | Add: linting_stack, structured_data, c_cpp_formatting |
| **Workflow schema** | Derive from the-rewrite-room frontmatter |
| **Reference agents** | Document: planning (swarm-task-planner, plan-validator), research (feature-researcher, codebase-analyzer), plugin_quality (plugin-assessor), documentation (plugin-docs-writer) |
| **Stack research** | Schema for references, research-curator population, research-context-agent integration |

---

## Cross-References

- **STATUS block**: the-rewrite-room, output-contracts.md — standardize across Layer 2 workflows
- **Fidelity rules**: summarizer — consider for summarization-type stacks
- **CLEAR+CoVe**: swarm-task-planner — task design standard for Layer 2 task plans
- **RT-ICA**: work-backlog-item, plugin-creator — pre-planning gate for Layer 2 additions
- **SAM artifacts**: plugin-creator, development-harness — .planning/harness/, STATE.md patterns
