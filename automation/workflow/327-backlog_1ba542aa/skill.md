---
last-updated: 2026-02-13
p0-count: 0
p1-count: 10
p2-count: 8
ideas-count: 10
---

# Backlog

Tracked features, ideas, and deferred work for grooming and future sessions.

---

## P0 - Must Have

_(Empty)_

---

## P1 - Should Have

### SAM: Error Recovery / Rollback Procedures

**Source**: Gap analysis of SAM framework
**Added**: 2026-02-01
**Description**: Define explicit procedure when a task fails irrecoverably. How to undo artifact changes? How to restore artifact plane to consistent state after failure?
**Research first**: How do GSD, BMAD-METHOD, AutoGPT, and traditional CI/CD handle rollback? What patterns exist for transactional artifact updates?
**Suggested location**: `methodology_development/stateless-software-engineering-framework.md` (new Appendix or Part 6 addition)

### SAM: Human Escalation Criteria

**Source**: Gap analysis of SAM framework
**Added**: 2026-02-01
**Description**: Define explicit triggers for escalating to human at each stage (not just Discovery). When should agents block and ask vs attempt repair vs fail?
**Research first**: How do GSD deviation rules work? How does BMAD-METHOD handle human checkpoints? What escalation patterns exist in agent frameworks?
**Suggested location**: `methodology_development/stateless-software-engineering-framework.md` (each Agent Specification section)

### SAM: Timeout/Stall Detection

**Source**: Gap analysis of SAM framework
**Added**: 2026-02-01
**Description**: Define mechanism to detect when an agent is stuck or has stalled. Include timeout thresholds per stage, health check patterns, and recovery actions.
**Research first**: How do orchestration frameworks (Temporal, Prefect, Airflow) handle task timeouts? What heartbeat patterns exist? How does Gas Town handle session recycling?
**Suggested location**: `methodology_development/stateless-software-engineering-framework.md` (Orchestrator section 3.8)

### SAM: Artifact Schema Validation

**Source**: Gap analysis of SAM framework
**Added**: 2026-02-01
**Description**: Define formal validation rules or JSON schemas for artifact formats. Currently only templates provided. Enable automated validation at stage boundaries.
**Research first**: How do GSD artifacts (STATE.md, ROADMAP.md) enforce structure? What validation approaches exist in BMAD-METHOD? JSON Schema vs YAML validation vs custom parsers?
**Suggested location**: `methodology_development/sam-artifact-schemas/` (new directory with schema files)

### SAM: Scope Creep Detection

**Source**: Gap analysis of SAM framework
**Added**: 2026-02-01
**Description**: Define mechanism to detect when execution diverges from plan. How does Forensic Review detect that the execution agent solved a different problem than planned?
**Research first**: How does GSD plan-checker detect deviation? What diff/comparison techniques exist? How do code review tools detect scope creep in PRs?
**Suggested location**: `methodology_development/stateless-software-engineering-framework.md` (section 3.6 Forensic Review)

### Resolve 48 pre-existing ty (Astral type checker) diagnostics

**Source**: `uv run ty check .` run during session 2026-02-11
**Added**: 2026-02-11
**Description**: ty v0.0.16 found 48 diagnostics across the codebase. Breakdown by category:
- **25 unresolved-import errors**: PEP 723 inline script dependencies (`httpx`, `anthropic`, `mcp`, `pandas`, `frontmatter`, `defusedxml`, etc.) not in the project venv. Also `plugin_validator`, `implementation_manager`, `file_metrics` — local module imports that ty can't resolve.
- **12 possibly-missing-attribute warnings**: `TempDoc | None` and `Path | None` unions not narrowed before attribute access in `find-temp-documentation.py` and `get-task-context.py`.
- **6 invalid-argument-type errors**: Functions receiving `T | None` when they expect `T` — same narrowing issue.
- **3 unsupported-operator / invalid-key errors**: Path `/` operator on union types, TypedDict subscript with runtime string key.
**Approach**:
1. Configure `ty.toml` or `pyproject.toml [tool.ty]` to exclude PEP 723 scripts (their deps are resolved at runtime, not in the venv)
2. Fix real type narrowing issues in `find-temp-documentation.py` and `get-task-context.py` (add `assert doc is not None` or `if doc is None: return` guards)
3. Fix `TypedDict` subscript issue with proper typing
4. Consider adding `ty check` to CI alongside basedpyright and mypy
**Files affected**: `.claude/utilities/find-temp-documentation.py`, `plugins/python3-development/skills/implementation-manager/scripts/get-task-context.py`, `plugins/plugin-creator/scripts/plugin-validator.py`, `plugins/summarizer/tests/test_file_metrics.py`

### Meta-Process Capture — Expert Panel Dataset Builder

**Source**: ARL expert panel process (sessions 2026-02-12 to 2026-02-13)
**Added**: 2026-02-13
**Description**: Document the multi-agent expert panel methodology as a reusable system for building datasets that inform skills and systems. The process — assign framework experts to repositories, ask structured questions, cross-examine, synthesize, map to requirements, validate — produced high-quality sourced findings. Capture this as a repeatable pattern.
**Key elements to document**:
- Expert assignment protocol (one agent per source repo, source-code-only evidence standard)
- Question group design (themed questions → cross-examination → synthesis)
- Cross-examination as adversarial validation (experts challenge each other's claims)
- Phased output (discussion → requirement mapping → synthesis → validation)
- Traceability chain (claim → expert citation → file:line evidence)
- Session continuity handling (state file enables cross-session resumption)
**Input artifacts**: `plugins/plugin-creator/skills/assessor/references/ARL/ARL-agent-instructions.md`, `qa-expert-panel.md`
**Suggested location**: New skill or methodology document — captures the meta-process, not the ARL content

### SAM Extension — Integrate ARL General Theory

**Source**: ARL expert panel Phase 3 output
**Added**: 2026-02-13
**Description**: Integrate the 7 universal principles from `synthesis-general-theory.md` into SAM methodology documents. These principles (structure over instruction, front-loading reduces gates, AI cannot self-evaluate, compression is architectural, iteration-aware state required, parallelism enables independent verification, failure paths need more compression) extend SAM's scope to cover autonomous refinement loops.
**Input artifacts**: `plugins/plugin-creator/skills/assessor/references/ARL/synthesis-general-theory.md`
**Target files**: `methodology_development/stateless-agent-methodology.md`, `methodology_development/stateless-software-engineering-framework.md`
**Dependencies**: None — general theory is framework-agnostic
**Related backlog items**: SAM gap items (error recovery, human escalation, scope creep detection) — the general theory findings inform several of these

### ARL Skill Development

**Source**: ARL expert panel Phase 3 output
**Added**: 2026-02-13
**Description**: Build the Autonomous Refinement Loop as a skill using `synthesis-arl-applicable.md` as its reference foundation. The ARL is a logical process (Assess → Plan → Implement → Review → Repeat) for autonomous skill refinement with R1-R10 requirement gates.
**Input artifacts**: `plugins/plugin-creator/skills/assessor/references/ARL/synthesis-arl-applicable.md`, `synthesis-general-theory.md`
**Scope**: Logical process design — what gates fire when, what each gate checks, success/failure criteria. NOT implementation artifacts (schemas, thresholds, pseudocode).
**Dependencies**: Benefits from SAM extension being done first (ARL would be a SAM-based skill)
**Suggested location**: New skill under `plugins/plugin-creator/skills/` or standalone plugin

### Extract claude-plugin-lint to standalone PyPI package

**Source**: Gap analysis - no existing Claude Code plugin linters exist
**Added**: 2026-02-01
**Description**: Extract and enhance `validate_frontmatter.py` into a standalone open-source project. First dedicated linter for Claude Code plugin frontmatter (SKILL.md, agents/*.md, commands/*.md). Official `claude plugin validate` only checks plugin.json structure.
**Features to include**:
- YAML frontmatter schema validation with Pydantic models
- Auto-fix capabilities (arrays → comma-separated, multiline → single-line)
- Token-based complexity metrics (tiktoken) instead of line counts
- Cross-reference validation (agent references non-existent skill)
- Marketplace readiness scoring
- Pre-commit hook integration
- CLI with `--fix` and `--report` modes
**Current source**: `plugins/plugin-creator/scripts/validate_frontmatter.py`
**Suggested repo name**: `claude-plugin-lint` or `cc-plugin-validator`

---

## P2 - Could Have

### SAM: Parser regex false positive on "## Task Summary Statistics"

**Source**: Migration proof-of-concept (2026-02-13)
**Added**: 2026-02-13
**Description**: The widened task header regex `^#{2,3}\s+Task:?\s+([A-Za-z0-9.]+)[:\s-]+(.+)$` in `implementation_manager.py` matches `## Task Summary Statistics` as task ID "Summary" with title "Statistics". The regex needs a negative lookahead or post-parse filter to exclude non-task sections. Observed when parsing `plan/tasks-1-plugin-linter.md`.
**File**: `plugins/python3-development/skills/implementation-manager/scripts/implementation_manager.py` line 645

### SAM: Replace validate-task-file.sh with Python validator

**Source**: Task format standardization plan (2026-02-13)
**Added**: 2026-02-13
**Description**: The bash validator at `plugins/plugin-creator/scripts/validate-task-file.sh` validates a different schema (`tasks-refactor-*.md`) and doesn't understand YAML frontmatter. Replace with Python validator that uses the shared `task_format.py` module.
**File**: `plugins/plugin-creator/scripts/validate-task-file.sh`

### SAM: Parallel Execution Details

**Source**: Gap analysis of SAM framework
**Added**: 2026-02-01
**Description**: Detail safe parallelization within SAM pipeline. When can tasks run in parallel? How to handle merge conflicts? Reference GSD wave execution pattern.
**Research first**: How does GSD wave execution work in detail? How do task orchestrators (Temporal, Prefect) handle parallel dependencies? What conflict resolution patterns exist?
**Suggested location**: `methodology_development/stateless-software-engineering-framework.md` (new section 2.4 or Appendix)

### SAM: Multi-Model Strategy

**Source**: Gap analysis of SAM framework
**Added**: 2026-02-01
**Description**: Define guidance for using different models for different agent types. E.g., cheaper/faster models for simple verification, stronger models for planning.
**Research first**: How do agent frameworks handle model selection? What cost/quality tradeoffs exist? How does Claude Code's haiku/sonnet/opus selection work?
**Suggested location**: `methodology_development/stateless-software-engineering-framework.md` (Implementation Roadmap or new Appendix)

### SAM: Audit Trail / Observability

**Source**: Gap analysis of SAM framework
**Added**: 2026-02-01
**Description**: Beyond artifacts, define logging/metrics/tracing guidance. How to diagnose pipeline issues? What telemetry to capture?
**Research first**: How do GSD and BMAD-METHOD handle logging? What observability patterns exist in agent frameworks? OpenTelemetry for LLM workflows?
**Suggested location**: `methodology_development/stateless-software-engineering-framework.md` (new Appendix I)

### SAM: Partial Success Handling

**Source**: Gap analysis of SAM framework
**Added**: 2026-02-01
**Description**: Define how to represent and handle partial task success. Task completes some DoD items but not all. How is this state represented in artifacts?
**Research first**: How do GSD checkpoints represent partial progress? How do CI/CD systems handle partial test passes? What state machine patterns exist?
**Suggested location**: `methodology_development/stateless-software-engineering-framework.md` (section 3.5 Execution Agent output)

### SAM: Context Size Management

**Source**: Gap analysis of SAM framework
**Added**: 2026-02-01
**Description**: Define explicit guidance for measuring and managing context size per agent. What's the target token budget? How to detect context pressure?
**Research first**: How do agent frameworks measure context usage? What token counting approaches exist? How does Claude Code handle context limits internally?
**Suggested location**: `methodology_development/stateless-software-engineering-framework.md` (section 2.1 or Appendix C)

### SAM: Conflicting Review Findings

**Source**: Gap analysis of SAM framework
**Added**: 2026-02-01
**Description**: Define protocol when forensic review and self-verification disagree. Which takes precedence? How to adjudicate conflicts?
**Research first**: How do code review systems handle conflicting reviewers? What adjudication patterns exist in multi-agent systems? How does GSD handle verification disagreements?
**Suggested location**: `methodology_development/stateless-software-engineering-framework.md` (section 3.6 Forensic Review)

---

## Ideas

### SAM: Cost/Token Management

**Source**: Gap analysis of SAM framework
**Added**: 2026-02-01
**Description**: Explore token budgets and cost controls per agent/stage. Track API costs. Set limits per task.
**Research first**: How do LLM cost management tools work (LangSmith, Helicone)? What budget enforcement patterns exist?

### SAM: Team Coordination Protocols

**Source**: Gap analysis of SAM framework
**Added**: 2026-02-01
**Description**: Explore how multiple humans interact with the pipeline concurrently. Locking? Ownership? Notifications?
**Research first**: How do collaborative editing systems handle concurrent users? What patterns exist in git-based workflows?

### SAM: External System Integration Patterns

**Source**: Gap analysis of SAM framework
**Added**: 2026-02-01
**Description**: Explore patterns for integrating with issue trackers (GitHub Issues, Jira), CI/CD pipelines, git hooks.
**Research first**: How does GSD integrate with external tools? What MCP servers exist for issue trackers? How do agent frameworks bridge to CI/CD?

### SAM: Migration Strategy Guide

**Source**: Gap analysis of SAM framework
**Added**: 2026-02-01
**Description**: Explore how to migrate existing projects to SAM. Incremental adoption? Parallel running? Artifact bootstrap?
**Research first**: How do organizations adopt new methodologies incrementally? What migration patterns exist for process changes?

### SAM: Training/Onboarding Materials

**Source**: Gap analysis of SAM framework
**Added**: 2026-02-01
**Description**: Explore creating materials for new team members learning the methodology. Tutorials, examples, quick-start guides.
**Research first**: What training materials exist for GSD, BMAD-METHOD? What makes effective agent methodology onboarding?

### SAM: Non-Code Workflow Guidance

**Source**: Gap analysis of SAM framework
**Added**: 2026-02-01
**Description**: Explore how SAM handles documentation-only tasks, configuration changes, infrastructure work. Adapt templates?
**Research first**: How do GSD and BMAD-METHOD handle non-code work? What artifact types exist beyond code?

### Carbonyl Browser Integration for Claude Code

**Source**: Session experimentation 2026-02-05
**Added**: 2026-02-05
**Description**: Research whether carbonyl (terminal Chromium browser) can work with Claude Code for reliable web content extraction. Carbonyl renders pages in terminal but needs a TTY.
**Research areas**:
- Can carbonyl run via tmux/screen/script to provide a pseudo-TTY?
- Could carbonyl be wrapped with a screenshot tool (e.g., termshot, asciinema) that passes images back to Claude?
- What's the minimal TTY setup needed for headless carbonyl operation?
- Compare with is-fast, lynx, w3m for text extraction capabilities
**Context**: WebFetch is unreliable (summarizing agents hallucinate), Playwright requires browser downloads that may be blocked. Carbonyl is self-contained but needs TTY.

### Validate is-fast for Web Content Extraction

**Source**: Session experimentation 2026-02-05
**Added**: 2026-02-05
**Description**: Test is-fast CLI tool on host with unrestricted network access.
**Validation steps**:
- Install: `curl --proto '=https' --tlsv1.2 -LsSf https://github.com/Magic-JD/is-fast/releases/latest/download/is-fast-installer.sh | sh`
- Test: `is-fast --direct https://code.claude.com/docs/en/skills --piped`
- Verify it extracts text content from JS-rendered pages
- Compare output quality with curl, lynx, w3m
- Test CSS selector filtering with `--selector`
**Blocked on 2026-02-05**: DNS resolution failed in restricted environment

### Validate agent-browser for Web Automation

**Source**: Session experimentation 2026-02-05
**Added**: 2026-02-05
**Description**: Test agent-browser (Playwright-based) on host with unrestricted network and Playwright browsers installed.
**Validation steps**:
- Install browsers: `npx playwright install`
- Test: `npx agent-browser open https://code.claude.com/docs/en/skills`
- Test: `npx agent-browser snapshot -i` (get element refs)
- Test: `npx agent-browser get text body` (extract page text)
- Verify snapshot/interact/re-snapshot workflow works
- Document prerequisites for skill to function
**Blocked on 2026-02-05**: Could not download Playwright browsers (DNS resolution failed, missing system libs)
**Skill location**: `.claude/skills/agent-browser/SKILL.md`

### Validate carbonyl Terminal Browser

**Source**: Session experimentation 2026-02-05
**Added**: 2026-02-05
**Description**: Test carbonyl on host with proper TTY and network access.
**Validation steps**:
- Test basic: `npx -y carbonyl --no-sandbox https://example.com`
- Test with tmux: `tmux new-session -d -s carbonyl 'npx -y carbonyl --no-sandbox https://example.com'`
- Test screenshot capture: Can we grab terminal output as image?
- Test text extraction: Can we pipe output or capture rendered text?
- Compare JS rendering quality with other tools
**Blocked on 2026-02-05**: Needs TTY (Inappropriate ioctl for device), DNS also blocked

---

## Completed

### Replace requests with httpx in all scripts

**Source**: CI/pre-commit inconsistency discovery (2026-02-05)
**Completed**: 2026-02-06
**Description**: Migrated `validate-glfm.py` from `requests` to `httpx`. Added TID251 ruff ban rule for `requests` imports. Removed `types-requests` from dev dependencies. `sync-gitlab-docs.py` was already using httpx.
**Location**: `plugins/gitlab-skill/skills/gitlab-skill/scripts/validate-glfm.py`, `pyproject.toml`

### Enhance swarm-task-planner with multi-source synthesis

**Source**: [external-pattern-integration-2026-02-01.md](.claude/external-pattern-integration-2026-02-01.md)
**Completed**: 2026-02-06
**Description**: Multi-source synthesis implemented via CLEAR+CoVe standard and Project Awareness section with investigation commands for searching documentation, assessing project structure, and identifying progress.
**Location**: `plugins/python3-development/agents/swarm-task-planner.md`

### Add context compliance checking

**Source**: [external-pattern-integration-2026-02-01.md](.claude/external-pattern-integration-2026-02-01.md)
**Completed**: 2026-02-06
**Description**: Complete plan-validator agent implements 8 validation dimensions including Requirement Coverage, Task Completeness, Dependency Correctness, Agent Capability Match, Input/Output Validity, Artifact Wiring, Testability, and Scope Sanity.
**Location**: `plugins/python3-development/agents/plan-validator.md`

### SAM: Artifact Versioning Strategy

**Source**: Gap analysis of SAM framework
**Completed**: 2026-02-06
**Description**: Implemented storage-agnostic semantic tokens using pattern `ARTIFACT:{TYPE}({SCOPE_OR_ID})` with disambiguators (CTX, PREREQ, EXEC, VERIFY). Provides both filesystem-backed and SQL-backed example implementations.
**Location**: `methodology_development/stateless-software-engineering-framework.md` (section 2.1.2)

### Create ecosystem-researcher agent

**Source**: [external-pattern-integration-2026-02-01.md](.claude/external-pattern-integration-2026-02-01.md)
**Completed**: 2026-02-05
**Description**: New agent for ecosystem/domain research before roadmap creation. Supports three modes - Ecosystem discovery, Feasibility assessment, Comparison analysis.
**Patterns from**: gsd-project-researcher.md (research modes)
**Location**: `plugins/python3-development/agents/ecosystem-researcher.md`

---

## Format Guide

```markdown
### Item title

**Source**: [link or description of where this came from]
**Added**: YYYY-MM-DD
**Description**: What needs to be done
**Research first**: (for SAM items) Questions to answer via competitive analysis before implementing
**Patterns from**: (optional) External source if from pattern integration
**Suggested location**: (optional) Where this should be implemented
```

### Research Resources

**Skill for research phase**: `/research-and-compare <methodology-name-or-url>`

- Produces structured comparison documents following SAM comparison template
- Includes overlap/divergence analysis, weakness discovery, implementation pairing
- Outputs to `methodology_development/.meta/v1_comparisons/`

**Existing SAM comparisons** (start here before running new research):

- [methodology_development/.meta/v1_comparisons/](../methodology_development/.meta/v1_comparisons/)
  - sam-vs-get-shit-done.md
  - sam-vs-bmad-method.md
  - sam-vs-gastown.md
  - sam-vs-taskmaster.md
  - sam-vs-octocode.md
  - sam-vs-superclaude.md
  - sam-vs-ralph-loop-orchestrator.md
  - sam-vs-cc-sessions.md
  - sam-vs-v-model.md
  - sam-infrastructure-layer.md

**Workflow for SAM gap items**:

1. Check existing comparisons for relevant findings
2. If more research needed: `/research-and-compare <framework>` for specific topics
3. Synthesize findings into SAM framework update
4. Mark backlog item complete

### P1: plugin-validator UX and coverage gaps

**Source**: Experimental validation of plugin-validator.py against all component types (2026-02-13)
**Added**: 2026-02-13
**Last groomed**: 2026-02-13
**File**: `./plugins/plugin-creator/scripts/plugin-validator.py` (2934+ lines)
**Tests**: `./plugins/plugin-creator/tests/` (12 test files, 93% pass rate per QA report)
**QA report**: `./plugins/plugin-creator/planning/plugin-validator-qa-report.md`

#### Sub-issue 1: UX — report counts validator results, not files

**Severity**: UX bug
**Lines**: 2928-2931 (result collection loop), 2698-2702 (summary display)

**Root cause**: Lines 2928-2931 iterate over `validators` (not files), appending `(path, result)` for each validator. When 1 file has 4 validators (FrontmatterValidator, NameFormatValidator, DescriptionValidator, NamespaceReferenceValidator), the `results` list has 4 entries all pointing to the same path. The report loop at line 2630 then prints "PASSED" 4 times for 1 file, and the summary at line 2702 shows "Total files: 4".

**Acceptance criteria**:
- Running validator on 1 file shows "Total files: 1"
- Each validator result is labeled with the validator name (e.g., "FrontmatterValidator: PASSED")
- Summary counts unique files, not validator invocations

#### Sub-issue 2: Commands receive skill-specific SK005 warning

**Severity**: False positive
**Lines**: 1884-1900 (SK005 check in DescriptionValidator), 2896-2901 (validator selection)

**Root cause**: `DescriptionValidator.validate()` (line 1802) has no file-type awareness. It applies SK004 ("description too short") and SK005 ("missing trigger phrases") to all file types. The validator selection at line 2896 applies DescriptionValidator to SKILL, AGENT, and COMMAND equally. Commands have a different frontmatter schema — they use `argument-hint`, `allowed-tools`, `agent` fields and do not need trigger phrases in their description.

**Acceptance criteria**:
- SK005 only fires on SKILL files (not COMMAND or AGENT)
- SK004 fires on SKILL and AGENT files (both need meaningful descriptions) but not COMMAND
- New error code series CM001+ for command-specific checks (e.g., validate `allowed-tools` format, `argument-hint` presence)
- DescriptionValidator receives file type context (via constructor parameter or path inspection)

#### Sub-issue 3: Hooks have zero validation

**Severity**: Missing feature
**Lines**: 148-165 (detect_file_type returns UNKNOWN for .js hooks and hooks.json)

**Root cause**: `FileType` enum (line 141-145) has no HOOK variant. `detect_file_type()` only matches SKILL.md, plugin.json, agents/*, commands/*. Hook files (`.js` in `hooks/` directories, `hooks.json` configs) fall through to UNKNOWN, which triggers the error at line 2920.

**Acceptance criteria**:
- `FileType` enum includes HOOK_SCRIPT and HOOK_CONFIG (or combined HOOK)
- `detect_file_type()` recognizes `.js` files in `hooks/` directories
- `detect_file_type()` recognizes `hooks.json` files
- New error code series HK001+ for hook validation (e.g., valid JSON in hooks.json, valid event types, valid hook matcher patterns)
- New `HookValidator` class following existing Validator protocol
- Reference: `/plugin-creator:claude-hooks-reference-2026` skill for hooks.json schema

#### Sub-issue 4: Dead code — nested skill resolution pattern

**Severity**: Code quality
**Lines**: 904-911 (in `_resolve_skill_reference` method of NamespaceReferenceValidator)

**Root cause**: Lines 904-911 search for `plugins/{plugin}/skills/*/{name}/SKILL.md` (nested category pattern). This pattern does not exist in the repository. Actual skill layout is `plugins/{plugin}/skills/{name}/SKILL.md` (flat). The direct check at lines 899-902 handles all real cases. The nested search at 904-911 is dead code. Additionally, lines 759-760 and 771-772 reference the nested pattern in error message strings.

**Note**: The original backlog item referenced lines 651-656 and 778-785 as dead code for `_discover_skills` and `_discover_invocable_skills`. Those methods do not exist in the current source. The actual dead code is in `_resolve_skill_reference` at lines 904-911.

**Acceptance criteria**:
- Remove lines 904-911 (nested skill resolution)
- Update error message strings at lines 759-760 and 771-772 to remove nested pattern references
- Verify no test depends on nested skill resolution (grep tests for `skills/*/` double-nested patterns)

#### Dependencies

- Sub-issues 1-4 are independent and can be implemented in parallel
- Sub-issue 3 (hooks) is the largest — requires new FileType, new Validator class, new error codes
- Sub-issue 2 (command SK005) requires deciding on a command-specific error code series

#### Related prior work

- QA report documents 15 pre-existing test failures (93% pass rate) — check for conflicts before modifying validators
- `/plugin-creator:claude-hooks-reference-2026` — hooks.json schema reference for sub-issue 3
- `/plugin-creator:claude-skills-overview-2026` — skill vs command schema differences for sub-issue 2

#### Approach

Use `/python3-development:add-new-feature` to extend `plugin-validator.py`. Delegate implementation to `@python-cli-architect`. Each sub-issue is a separate task.
