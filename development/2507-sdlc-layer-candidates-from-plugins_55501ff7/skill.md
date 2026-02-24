# SDLC Layer Separation Architecture — Plugin Candidates

**Generated**: 2026-02-23
**Source**: Full read of all plugins under `plugins/` (CLAUDE.md, plugin.json, SKILL.md, agents, references, workflows)
**Plan Context**: Layer 0 (SDLC-Agnostic), Layer 1 (Language-Specific), Layer 2 (Stack/Goal-Specific), ARL Meta-Layer

---

## Layer 0 (SDLC-Agnostic)

### development-harness (plugin)
- **Type**: methodology
- **Path**: `plugins/development-harness/CLAUDE.md`
- **Proposed Layer**: L0
- **Relevance**: Core SAM 7-stage pipeline implementation; language-agnostic process harness with ARL human touchpoints and Voltron-style composition. Defines artifact conventions, stage gating, and human escalation model.
- **Key content**: S1 Discovery → S2 Planning+RT-ICA → S3 Context Integration → S4 Task Decomposition → S5 Execution → S6 Forensic Review → S7 Final Verification; ARTIFACT:{TYPE}({ID}) token pattern; `.planning/harness/` layout; harness owns process, language plugins own specialists.

### development-harness/default-development-flow
- **Type**: workflow
- **Path**: `plugins/development-harness/skills/development-harness/references/default-development-flow.md`
- **Proposed Layer**: L0
- **Relevance**: Canonical SAM pipeline with ARL touchpoint gates; stage descriptions, artifact types, loop limits (3 NEEDS_WORK, 2 NOT_CERTIFIED).
- **Key content**: Mermaid flowchart of pipeline; Gate 1 (S1→S2) for unbound constraints; Gate 2 (S4→S5) for high complexity; flow override rules for language plugins.

### development-harness/human-touchpoint-model
- **Type**: reference
- **Path**: `plugins/development-harness/skills/development-harness/references/human-touchpoint-model.md`
- **Proposed Layer**: L0
- **Relevance**: ARL-derived escalation decision model; when to escalate vs proceed autonomously.
- **Key content**: Bound/Unbound/Mixed constraint types; risk assessment (reversible, local scope, existing patterns); domain knowledge assessment; pre-scheduled gates; dynamic escalation (NEEDS_WORK limit, NOT_CERTIFIED limit).

### development-harness/artifact-conventions
- **Type**: reference
- **Path**: `plugins/development-harness/skills/development-harness/references/artifact-conventions.md`
- **Proposed Layer**: L0
- **Relevance**: SAM artifact naming, file layout, cross-referencing for stateless handoff.
- **Key content**: ARTIFACT:{TYPE}({SCOPE_OR_ID}) pattern; required sections per artifact type (DISCOVERY, PLAN, CONTEXT, TASK, EXECUTION, REVIEW, VERIFICATION); `.planning/harness/` coexistence with other tools.

### development-harness/role-resolution-protocol
- **Type**: reference
- **Path**: `plugins/development-harness/skills/development-harness/references/role-resolution-protocol.md`
- **Proposed Layer**: L0
- **Relevance**: How abstract roles (architect, test-designer, code-reviewer, design-spec, linting) resolve to concrete agents at runtime.
- **Key content**: Project detection markers (pyproject.toml, package.json, Cargo.toml); manifest search order; fallback to general-purpose; quality gate loading.

### development-harness/language-manifest-schema
- **Type**: reference
- **Path**: `plugins/development-harness/skills/development-harness/references/language-manifest-schema.md`
- **Proposed Layer**: L0
- **Relevance**: Schema for language plugins to declare Role Fulfillment, Quality Gates, Project Detection, optional Flow Override.
- **Key content**: Role Fulfillment format; Quality Gates (format, lint, typecheck, test, standards); Project Detection (markers, source-patterns, test-patterns); Process Flow Override rules.

### development-harness/workflow stages (Discovery, Planning, Context Integration, Task Decomposition, Execution, Forensic Review, Final Verification)
- **Type**: skill
- **Path**: `plugins/development-harness/skills/workflows/*/SKILL.md`
- **Proposed Layer**: L0
- **Relevance**: Seven SAM stage implementations; each produces named artifact per stage.
- **Key content**: Discovery (WHO/WHAT/WHEN/WHY, no HOW); Planning (RT-ICA pre-pass); Task Decomposition (task files with acceptance criteria); Execution (delegate to specialists); Forensic Review (COMPLETE/NEEDS_WORK); Final Verification (CERTIFIED/NOT_CERTIFIED).

### development-harness/planner-rt-ica
- **Type**: skill
- **Path**: `plugins/development-harness/skills/planner-rt-ica/SKILL.md` (via python3-development)
- **Proposed Layer**: L0
- **Relevance**: Planning-phase Input Completeness Analysis; prevents invented requirements while allowing dependency-first plan.
- **Key content**: APPROVED-FOR-PLANNING / APPROVED-WITH-GAPS / BLOCKED-FOR-PLANNING; Missing Inputs by dependency; Required Unblock Actions; Planning Annotations; no invention rule.

### development-harness/validation-protocol, clear-cove-task-design, generate-task
- **Type**: skill
- **Path**: `plugins/development-harness/skills/`
- **Proposed Layer**: L0
- **Relevance**: Planning tools: validation patterns, task design methodology, task file generation.
- **Key content**: Validation checklists; CoVe-style task design; task file format with inputs, acceptance criteria, agent assignment.

### verification-gate
- **Type**: skill
- **Path**: `plugins/verification-gate/skills/verification-gate/SKILL.md`
- **Proposed Layer**: L0
- **Relevance**: Mandatory pre-action verification checkpoints; prevents pattern-matching from overriding explicit reasoning. Blocks execution when hypothesis unverified or action targets different system.
- **Key content**: 4 checkpoints (Hypothesis Stated, Hypothesis Verified, Hypothesis-Action Alignment, Pattern-Matching Detection); VERIFICATION COMPLETE / EXECUTION BLOCKED format; CoVe integration; research-foundations.md, failure-patterns.md references.

### orchestrator-discipline
- **Type**: practice
- **Path**: `plugins/orchestrator-discipline/CLAUDE.md`, `plugins/orchestrator-discipline/rules/CLAUDE.md`, `plugins/orchestrator-discipline/skills/orchestrator-discipline/SKILL.md`
- **Proposed Layer**: L0
- **Relevance**: Context window discipline for orchestrator role; read constraints, delegation constraints, anti-patterns.
- **Key content**: PERMITTED vs NEVER reads (source files, config, test files, diagnostic output, agent transcripts); Investigation Escalation Anti-Pattern (4-step sequence); Agent Output Polling Anti-Pattern; Diagnostic Commands delegation; Falsifiable test ("Will I Edit/Write this file this turn?"); PreToolUse hooks.

### agent-orchestration
- **Type**: skill
- **Path**: `plugins/agent-orchestration/skills/agent-orchestration/SKILL.md`
- **Proposed Layer**: L0
- **Relevance**: Scientific delegation framework; WHERE-WHAT-WHY context patterns; preserve agent autonomy.
- **Key content**: Observations + Success Criteria + Available Resources - Assumptions - Prescriptions; Delegation Template (OBSERVATIONS, DEFINITION OF SUCCESS, CONTEXT, YOUR TASK); ECOSYSTEM CONTEXT rules (session-specific only); Pre-Gathering Anti-Pattern; Holistic vs Micromanaged delegation; Pattern Expansion (single instance → systemic fix).

### agent-orchestration/how-to-delegate
- **Type**: skill
- **Path**: `plugins/agent-orchestration/skills/how-to-delegate/SKILL.md`
- **Proposed Layer**: L0
- **Relevance**: 10-step delegation worksheet; pre-flight verification before Task tool.
- **Key content**: Step 1 Load Framework; Step 2 Identify Task Type; Step 3 Gather Observations (no pre-gathering); Step 4 Define Success; Step 5 World-Building Context; Step 6 Available Resources; Step 7 Select Agent; Step 8 Pre-Flight Verification; Step 9 Construct Task Prompt; Step 10 Delegate.

### hallucination-detector
- **Type**: practice
- **Path**: `plugins/hallucination-detector/commands/hallucination-audit.md`, `plugins/hallucination-detector/hooks.json`
- **Proposed Layer**: L0
- **Relevance**: Stop-hook hallucination detector; blocks stopping when speculation-as-diagnosis, invented causality, pseudo-quantification, completeness overclaims detected.
- **Key content**: 5 triggers (Speculation Language, Causality Without Evidence, Pseudo-Quantification, Completeness Claims, Delegation Micromanagement); Pass/Fail; rewrite requirements; PreStop hook via node script.

### the-rewrite-room/output-contracts
- **Type**: reference
- **Path**: `plugins/the-rewrite-room/skills/the-rewrite-room/references/output-contracts.md`
- **Proposed Layer**: L0
- **Relevance**: STATUS block format for workflow termination; DONE/BLOCKED/FAILED; artifact listing; validation results.
- **Key content**: status-block-v1 (STATUS, SUMMARY, ARTIFACTS, VALIDATION, NOTES); audit-block-v1 (evidence requirements); summary-block-v1; citation-check constraint.

### the-rewrite-room/fidelity-rules
- **Type**: reference
- **Path**: `plugins/the-rewrite-room/skills/the-rewrite-room/references/fidelity-rules.md`
- **Proposed Layer**: L0
- **Relevance**: Anti-hallucination rules for summarization; read before summarizing, extract before abstracting, preserve counts.
- **Key content**: Rule 1 Read Before Summarizing; Rule 2 Extract Before Abstracting; Rule 3 Preserve Counts; Rule 4 Distinguish Absence from Nonexistence; Rule 5 No Lossy Re-Summarization; Rule 6 State Confidence; Rule 7 Structured Output Always.

### plugin-creator/arl
- **Type**: skill
- **Path**: `plugins/plugin-creator/skills/arl/SKILL.md`
- **Proposed Layer**: L0
- **Relevance**: Autonomous Refinement Loop knowledge reference; HOOTL prerequisites; 10 gates (R1-R10); when human gates can be replaced.
- **Key content**: HOOTL concept; Three Layers (Research Body, Execution Model, Observation); ARL-SAM-kaizen relationship triangle; 10 Gates (Information Completeness, Loop Detection, Validity Filtering, Plan Quality, Purpose Anchor, Content-Loss Detection, Convergence Tracking, Proportionality Check, Downstream Impact, Split Justification); Decision Tree (4 conditions); Universal Principles; Scope-Feasibility Matrix.

---

## Layer 1 (Language-Specific)

### development-harness/language-manifest-template
- **Type**: reference
- **Path**: `plugins/development-harness/templates/language-manifest-template.md`
- **Proposed Layer**: L1
- **Relevance**: Template for language plugins to declare Role Fulfillment, Quality Gates, Project Detection.
- **Key content**: Skeleton for Python, TypeScript, Rust, Go manifests; abstract role mapping.

### python3-development
- **Type**: skill
- **Path**: `plugins/python3-development/skills/python3-development/SKILL.md`
- **Proposed Layer**: L1
- **Relevance**: Python 3.11+ language manifest; orchestration guide, agent selection, quality gates, project structure.
- **Key content**: ROLE_TYPE identification; Linting Discovery Protocol; Format-First; Type Checker Discovery; Quality Gates (ruff, mypy, pytest); Standard Project Structure (packages/, tests/); Agent orchestration (python-cli-architect, python-pytest-architect, python-code-reviewer); Pre-Delegation Checklist; python-development-orchestration.md reference.

### python3-development/references
- **Type**: reference
- **Path**: `plugins/python3-development/skills/python3-development/references/python-development-orchestration.md`, `user-project-conventions.md`, `tool-library-registry.md`
- **Proposed Layer**: L1
- **Relevance**: Python-specific workflow patterns, conventions, tool catalog.
- **Key content**: TDD, Feature Addition, Code Review, Refactoring, Debugging workflows; agent chaining; multi-agent patterns; pyproject.toml template variables.

### perl-development
- **Type**: skill
- **Path**: `plugins/perl-development/` (skills: perl-testing, perl-cpan-ecosystem, perl-development, perl-environment-setup, perl-lint, perl-validate)
- **Proposed Layer**: L1
- **Relevance**: Perl 5.30+ language manifest; scripting, CPAN, testing, linting.
- **Key content**: perl-script-developer, perl-script-auditor, perl-cli-architect agents; modular skills for ecosystem, environment, validation.

### bash-development
- **Type**: skill
- **Path**: `plugins/bash-development/` (skills: bash-51-features, bash-52-features, bash-53-features, bash-development, bash-lint, bash-logging, bash-portability, bash-testing)
- **Proposed Layer**: L1
- **Relevance**: Bash 5.1+ and POSIX shell language manifest; version-specific features, portability, testing.
- **Key content**: bash-script-developer, bash-script-auditor agents; bash-logging, bash-portability, bash-testing skills.

---

## Layer 2 (Stack/Goal-Specific)

### the-rewrite-room
- **Type**: skill
- **Path**: `plugins/the-rewrite-room/skills/the-rewrite-room/SKILL.md`
- **Proposed Layer**: L2
- **Relevance**: Documentation workflow router; maps task intent to canonical workflows (drift-audit, documentation-sync, authoring, prompt-optimization, summarization, formatting-validation, research-utilities).
- **Key content**: Workflow taxonomy; registry/workflows.yaml; adapters for legacy components; output contract (STATUS block); router.py classify scripts.

### the-rewrite-room/workflows
- **Type**: workflow
- **Path**: `plugins/the-rewrite-room/skills/the-rewrite-room/workflows/*.md`
- **Proposed Layer**: L2
- **Relevance**: Canonical workflows for each documentation task type.
- **Key content**: drift-audit, documentation-sync, authoring, prompt-optimization, summarization, formatting-validation, research-utilities; adapters (summarizer, doc-drift-auditor, contextual-ai-documentation-optimizer, add-doc-updater).

### summarizer
- **Type**: skill
- **Path**: `plugins/summarizer/skills/summarizer/SKILL.md`
- **Proposed Layer**: L2
- **Relevance**: Summarization routing; format selection (structured, bullets, tldr, json, table, outline); decision tree for file/URL/image/multi-source.
- **Key content**: Fidelity Rules (mandatory); Team Coordination for 3+ sources; agent-result-relay; structured output with YAML frontmatter; anti-patterns.

### holistic-linting
- **Type**: skill
- **Path**: `plugins/holistic-linting/`
- **Proposed Layer**: L2
- **Relevance**: Code quality orchestration; format → lint → resolve via concurrent agents; linting-root-cause-resolver.
- **Key content**: holistic-linting-orchestrator, holistic-linting-resolver; post-linting-architecture-reviewer; ruff, mypy, bandit rules.

### gitlab-skill
- **Type**: skill
- **Path**: `plugins/gitlab-skill/skills/gitlab-skill/`
- **Proposed Layer**: L2
- **Relevance**: GitLab CI/CD, GLFM, gitlab-ci-local; pipeline debugging, GLFM syntax.
- **Key content**: .gitlab-ci.yml configuration; gitlab-ci-local-guide.md; validate_glfm.py; CI/CD optimization.

### fastmcp-creator
- **Type**: skill
- **Path**: `plugins/fastmcp-creator/skills/fastmcp-creator/`
- **Proposed Layer**: L2
- **Relevance**: MCP server creation; FastMCP 3.x framework; generic MCP protocol.
- **Key content**: references (community-practices, typescript-mcp-server, development-guidelines, mcp-best-practices, evaluation-guide, example-projects, prompts-and-templates).

### dasel
- **Type**: skill
- **Path**: `plugins/dasel/`
- **Proposed Layer**: L2
- **Relevance**: Structured data query/transform (JSON, YAML, TOML, XML, CSV, HCL, INI) via dasel v3.
- **Key content**: data-explorer, data-analyst, dasel-guide agents; data-exploration, data-transformation, setup skills; enterprise domain skills (hibernate, maven, spring, tomcat, installanywhere).

### clang-format
- **Type**: skill
- **Path**: `plugins/clang-format/skills/clang-format/`
- **Proposed Layer**: L2
- **Relevance**: C/C++ formatting configuration; .clang-format; brace styles, indentation, alignment.
- **Key content**: clang-format configuration patterns; preserve existing style; codify dominant conventions.

### xdg-base-directory
- **Type**: skill
- **Path**: `plugins/xdg-base-directory/skills/xdg-base-directory/`
- **Proposed Layer**: L2
- **Relevance**: Cross-platform config/data/cache/state file storage; XDG Base Directory spec; platformdirs.
- **Key content**: When to use; ~/.appname vs XDG; platformdirs implementation.

### plugin-creator
- **Type**: skill
- **Path**: `plugins/plugin-creator/CLAUDE.md`, `plugins/plugin-creator/skills/`
- **Proposed Layer**: L2
- **Relevance**: Plugin development toolkit; create/refactor/validate plugins, agents, skills; agent-creator, skill-creator, refactor workflows.
- **Key content**: plugin_validator.py (token-based complexity, frontmatter, links); create_plugin.py; agent-creator (scope, templates); refactor-plugin, refactor-skill; claude-skills-overview-2026, claude-plugins-reference-2026, claude-hooks-reference-2026; hooks-core-reference, hooks-patterns; add-doc-updater; auto_sync_manifests.py (version bumping).

### plugin-creator/agent-creator
- **Type**: skill
- **Path**: `plugins/plugin-creator/skills/agent-creator/SKILL.md`
- **Proposed Layer**: L2
- **Relevance**: Create agents from scratch or templates; scope (project/user/plugin); subagent-contract for role-based agents.
- **Key content**: agent-schema.md; discovery phase; template selection; scope determination; plugin.json update; validation.

### development-harness/agents (generic)
- **Type**: agent
- **Path**: `plugins/development-harness/agents/*.md`
- **Proposed Layer**: L2
- **Relevance**: Generic SAM workflow agents (swarm-task-planner, plan-validator, feature-researcher, codebase-analyzer, feature-verifier, integration-checker, context-gathering, context-refinement, doc-drift-auditor, ecosystem-researcher, service-docs-maintainer).
- **Key content**: subagent-contract; stateless artifact outputs; role-specific workflows.

---

## ARL Meta-Layer

### agentskill-kaizen
- **Type**: skill
- **Path**: `plugins/agentskill-kaizen/`
- **Proposed Layer**: ARL
- **Relevance**: Observe layer — session transcript analysis; identify inefficiencies, anti-patterns, repeated mistakes, missing tooling opportunities, user frustration signals.
- **Key content**: transcript-analysis skill (10 dimensions: Tool Misuse, Repeated Errors, User Frustration, Missing Tooling, Subagent Delegation, Shortest Path, Red Herring, System Interruptions, Missing Hooks, DuckDB SQL); kaizen-improvement skill (hook generation, agent patches, skill patches, CLAUDE.md updates, script automation); PM4Py process mining; kaizen-analysis MCP; improvement-generator, transcript-analyst agents.

### agentskill-kaizen/transcript-analysis
- **Type**: skill
- **Path**: `plugins/agentskill-kaizen/skills/transcript-analysis/SKILL.md`
- **Proposed Layer**: ARL (Observe)
- **Relevance**: JSONL schema, DuckDB query patterns, 10 analysis dimensions, process mining methodology.
- **Key content**: ~/.claude/projects/ transcript location; extract_tool_sequences, discover_process_model, check_conformance, find_frequent_patterns, detect_frustration_signals, cluster_sessions; output to .planning/kaizen/analysis-DATE.md.

### agentskill-kaizen/kaizen-improvement
- **Type**: skill
- **Path**: `plugins/agentskill-kaizen/skills/kaizen-improvement/SKILL.md`
- **Proposed Layer**: ARL (Improve)
- **Relevance**: Transform findings into actionable improvements; hooks, agent prompts, skill patches, CLAUDE.md, scripts.
- **Key content**: Hook generation (PreToolUse, SubagentStart, Stop); agent prompt refinement via subagent-refactorer; skill patches via skill-creator; output to .planning/kaizen/improvements/; --install flag for hooks.

### hallucination-detector
- **Type**: practice
- **Path**: `plugins/hallucination-detector/`
- **Proposed Layer**: ARL (Identify)
- **Relevance**: Identify deviations/hallucinations — speculation-as-diagnosis, invented causality, pseudo-quantification, completeness overclaims.
- **Key content**: Stop hook blocks; audit triggers; evidence-first rewrite requirement.

### plugin-creator/arl
- **Type**: skill
- **Path**: `plugins/plugin-creator/skills/arl/SKILL.md`
- **Proposed Layer**: ARL
- **Relevance**: ARL theory — failure categories, prerequisites, conditions; when human gates replaceable; scope-feasibility matrix.
- **Key content**: Layer 3 Observation (agentskill-kaizen as post-hoc implementation); AI user representatives; asynchronous feedback queue; question-to-action-item conversion.

### summarizer/session-historian
- **Type**: concept
- **Path**: (session-historian is in .claude/skills/ but summarizer supports session recall)
- **Proposed Layer**: ARL (Observe)
- **Relevance**: Session recall via transcript analysis; prior decisions, experiments, outcomes.
- **Key content**: DuckDB index; JSONL transcripts; session summaries.

---

## Summary by Plugin

| Plugin | Key L0 | Key L1 | Key L2 | Key ARL |
|--------|--------|--------|--------|---------|
| development-harness | SAM pipeline, artifact conventions, human touchpoints, role resolution, workflow stages | — | Generic agents | — |
| orchestrator-discipline | Read constraints, delegation anti-patterns, PreToolUse hooks | — | — | — |
| verification-gate | 4-checkpoint verification, hypothesis-action alignment | — | — | — |
| agent-orchestration | Delegation template, scientific method, pre-gathering anti-pattern | — | — | — |
| hallucination-detector | Stop-hook speculation detection | — | — | Identify |
| plugin-creator | ARL skill | — | Plugin creation, agent-creator, refactor | ARL theory |
| the-rewrite-room | Output contracts, fidelity rules | — | Workflow router, adapters | — |
| summarizer | Fidelity rules | — | Format routing, decision tree | — |
| python3-development | planner-rt-ica | Language manifest, orchestration guide | Quality gates, agents | — |
| perl-development | — | Language manifest | Perl agents, skills | — |
| bash-development | — | Language manifest | Bash agents, skills | — |
| agentskill-kaizen | — | — | — | Observe, Identify, Improve |
| holistic-linting | — | — | Lint orchestration | — |
| gitlab-skill | — | — | GitLab CI, GLFM | — |
| fastmcp-creator | — | — | MCP server creation | — |
| dasel | — | — | Data query/transform | — |
| clang-format | — | — | C/C++ formatting | — |
| xdg-base-directory | — | — | Config storage | — |

---

## Plugins Not Listed (Minimal SDLC Relevance)

- **brainstorming-skill** — Creative ideation; not SDLC pipeline
- **commitlint** — Commit message linting; L2 tool
- **conventional-commits** — Commit format; L2 convention
- **litellm** — LLM API proxy; infrastructure
- **llamafile** — Local model runner; infrastructure
- **prompt-optimization-claude-45** — Prompt optimization; L2 content
- **uv** — Python package manager; L2 tool
