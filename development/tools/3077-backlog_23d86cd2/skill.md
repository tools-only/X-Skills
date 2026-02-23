---
last-updated: 2026-02-22
last-completed: 2026-02-22
p0-count: 0
p1-count: 22
p2-count: 26
ideas-count: 11
---

# Backlog

Tracked features, ideas, and deferred work for grooming and future sessions.

---

## P0 - Must Have

_(Empty)_

---

## P1 - Should Have

### gitlab-skill: Remove hardcoded corporate URL

**Source**: Plugin code review session 2026-02-21
**Added**: 2026-02-21
**Completed**: 2026-02-22
**Priority**: P1
**Status**: DONE — replaced `https://sourcery.assaabloy.net` with generic placeholder `https://gitlab.example.com` in validate_glfm.py (default arg + help text) and gitlab-ci-local-guide.md (example URL)
**Description**: `validate_glfm.py` lines 152-153 hardcode `https://sourcery.assaabloy.net` as the default GitLab instance URL. `gitlab-ci-local-guide.md` line 51 also references this URL. This leaks a corporate internal URL into a public repository. Replace with a generic placeholder (e.g., `https://gitlab.example.com`) or make the URL a required argument with no default.
**Files**:
- `plugins/gitlab-skill/skills/gitlab-skill/scripts/validate_glfm.py` (lines 152-153)
- `plugins/gitlab-skill/skills/gitlab-skill/references/gitlab-ci-local-guide.md` (line 51)

### bash-development: Fix bash-53-features inaccuracies and task_output bug

**Source**: Plugin code review session 2026-02-21
**Added**: 2026-02-21
**Completed**: 2026-02-22
**Priority**: P1
**Status**: DONE — Fixed GLOBSORT syntax (`:asc`/`:desc` → `+`/`-` prefix, `date` → `mtime`, added missing specifiers), added required space in `${ command; }` and `${| command; }` examples, added REPLY-is-local note, corrected C23 claim to "C standard conformance improvements", added interface disclaimer for `kv`/`strptime`. Fixed `task_output()` bug in `log_functions.sh` — `${task_output}` → `${raw_task_output}` on lines 1410/1412.
**Description**: Two issues:

1. **bash-53-features/SKILL.md inaccuracies** (features exist but details wrong):
   - GLOBSORT: `:asc`/`:desc` suffixes fabricated — actual syntax is `+`/`-` prefix; `date` not a valid specifier (should be `mtime`); missing specifiers: `blocks`, `atime`, `ctime`, `numeric`, `nosort`
   - `${ command; }` examples: missing required space after `{` (bash.1 requires space/tab/newline/`|` after `{`)
   - `${| command; }` REPLY examples: misleading — REPLY is local within substitution, restored after completion
   - C23 claim overstated: build minimum is C90, not C23; C23 changes are about conformance
   - `kv`/`strptime` usage examples: unverifiable from official sources (existence confirmed, interface speculative)
2. **log_functions.sh bug** (line 1401): `task_output()` function references `${task_output}` variable on lines 1410/1412, but only `raw_task_output` is assigned (line 1406). Output will be empty.

**Citations**:
- Bash CHANGES: <https://tiswww.case.edu/php/chet/bash/CHANGES> lines 850-860 (accessed 2026-02-21)
- Bash NEWS: <https://ftp.gnu.org/gnu/bash/bash-5.3.tar.gz> extracted NEWS lines 48-57 (accessed 2026-02-21)
- Bash manpage: `bash-5.3/doc/bash.1` lines 2525-2566 (GLOBSORT), 4146+ (command substitution) (accessed 2026-02-21)
**Files**:
- `plugins/bash-development/skills/bash-53-features/SKILL.md` (GLOBSORT, examples, C23)
- `plugins/bash-development/skills/bash-logging/scripts/log_functions.sh` (lines 1401, 1410, 1412)

### ~~commitlint: Verify --last flag and exit codes against primary sources~~ RESOLVED

**Source**: Plugin code review session 2026-02-21
**Added**: 2026-02-21
**Resolved**: 2026-02-21
**Priority**: P1
**Status**: FACT-CHECKED 2026-02-21 — `--last` flag VERIFIED across 4 independent sources. Flag exists in source code (`cli.ts` with alias `-l`), official CLI reference docs, commitlint help output, and raw documentation. Exit codes verified against `ExitCode` enum in `cli-error.ts`. Original review agent claim that `--last` was fabricated was wrong.
**Citations**:
- commitlint source: `@commitlint/cli/src/cli.ts` (accessed 2026-02-21)
- commitlint docs: <https://commitlint.js.org/reference/cli.html> (accessed 2026-02-21)
- commitlint source: `@commitlint/cli/src/cli-error.ts` `ExitCode` enum (accessed 2026-02-21)

### clang-format: Fix broken YAML frontmatter

**Source**: Plugin code review session 2026-02-21
**Added**: 2026-02-21
**Completed**: 2026-02-22
**Priority**: P1
**Status**: DONE — Fixed `description:"Configure..."` → `description: "Configure..."` in `plugins/clang-format/skills/clang-format/SKILL.md` line 3.
**Description**: SKILL.md line 3 has `description:"Configure clang-format..."` — missing the required space after the `description:` key. This causes YAML parsing failures. The frontmatter should be `description: "Configure clang-format..."`.
**Files**:
- `plugins/clang-format/skills/clang-format/SKILL.md` (line 3)

### agent-orchestration: Remove phantom /is-it-done command references

**Source**: Plugin code review session 2026-02-21
**Added**: 2026-02-21
**Completed**: 2026-02-22
**Priority**: P1
**Status**: DONE — Replaced all 9 actionable `/is-it-done` slash command calls in `SKILL.md` and `how-to-delegate/SKILL.md` with `/am-i-complete` (the actual existing command). Source attribution references in `post-completion-validation-protocol.md` and `synthesis-improvements-from-research.md` left as historical metadata. Removed Clavix references from `clear-framework.md` — replaced with generic imperative descriptions of the CLEAR framework operations.

### perl-development: Fix shell injection vulnerability in example template

**Source**: Plugin code review session 2026-02-21
**Added**: 2026-02-21
**Completed**: 2026-02-22
**Priority**: P1
**Status**: DONE — Added `#!/usr/bin/env perl` shebang. Replaced undefined `App::Logger` import with inline logging functions (`print_success`, `print_error`, `print_info`, `print_debug`, `print_warning`). Added `use File::Spec`, `use IPC::Open3`, `use Symbol 'gensym'` for safe system calls. Fixed `command_exists` to use `File::Spec->path()` instead of backtick `command -v` shell call. Fixed `run_command` to accept list args and use `IPC::Open3::open3` instead of backtick string injection. Fixed `system("stty ...")` calls in `query_terminal_safely` to use list form.
**Files**:
- `plugins/perl-development/skills/perl-development/references/perl_example_file.pl`

### hallucination-detector: Fix dead backtick evidence marker and multi-occurrence "because" bug

**Source**: Plugin code review session 2026-02-21
**Added**: 2026-02-21
**Completed**: 2026-02-22
**Priority**: P1
**Status**: DONE — Fixed both bugs in `hallucination-audit-stop.js`: (1) Removed dead backtick regex from `EVIDENCE_MARKERS` (was never matched since `stripLowSignalRegions` removes backtick spans before scanning); replaced with a `BACKTICK_RE` check against the original pre-strip text via new `rawText` param in `hasEvidenceNearby`. (2) Changed `because` loop from `lower.indexOf()` (first occurrence only) to a `while` loop scanning all occurrences — each unsuppressed `because` is now flagged independently.
**Files**:
- `plugins/hallucination-detector/scripts/hallucination-audit-stop.js`

### python3-development: Fix 3 malformed frontmatter files

**Source**: Plugin code review session 2026-02-21
**Added**: 2026-02-21
**Completed**: 2026-02-22
**Priority**: P1
**Status**: DONE — Fixed 2 malformed frontmatter files (backlog said 3, validator found 2): `skills/development/add-new-feature/SKILL.md` and `skills/development/complete-implementation/SKILL.md` — both had `description:"..."` missing the space after the colon. Also unquoted the description in `complete-implementation` (was double-escaped `"\"..."\""`). No ghost agent reference found after inspection.

### the-rewrite-room: Fix nonexistent script reference and 6 broken links

**Source**: Plugin code review session 2026-02-21
**Added**: 2026-02-21
**Completed**: 2026-02-22
**Priority**: P1
**Status**: DONE — Fixed `validators.yaml` line 12: changed `validate_frontmatter.py` (no longer exists) to `plugin_validator.py`. Fixed 6 broken links in `registry-guide.md` (lines 9, 10, 11, 35, 57, 83) — updated `./workflows.yaml`, `./validators.yaml`, `./routing-rules.yaml` to `../registry/workflows.yaml`, `../registry/validators.yaml`, `../registry/routing-rules.yaml`. Added `--json` flag to `file_metrics.py` (`count` and `scan` commands) to support machine-readable output as referenced in `research-utilities.md`. Also fixed `validators.yaml` line 65 false-positive link pattern.
**Files**:
- `plugins/the-rewrite-room/skills/the-rewrite-room/registry/validators.yaml`
- `plugins/the-rewrite-room/skills/the-rewrite-room/references/registry-guide.md`
- `plugins/the-rewrite-room/skills/the-rewrite-room/workflows/research-utilities.md`
- `plugins/the-rewrite-room/skills/the-rewrite-room/scripts/file_metrics.py`

### fastmcp-creator: Add citations for 1200+ lines of FastMCP 3.x API documentation

**Source**: Plugin code review session 2026-02-21
**Added**: 2026-02-21
**Completed**: 2026-02-22
**Priority**: P1
**Status**: DONE — Added `## Sources` sections to all 7 reference files in `plugins/fastmcp-creator/skills/fastmcp-creator/references/`: `development-guidelines.md` (FastMCP GitHub, docs, PyPI, launch post, MCP spec, verified auth source files), `community-practices.md` (FastMCP GitHub, docs, PyPI, MCP spec), `evaluation-guide.md` (FastMCP GitHub, docs, MCP spec), `example-projects.md` (FastMCP GitHub, docs, MCP server registry, TypeScript SDK), `mcp-best-practices.md` (FastMCP GitHub, docs, MCP spec, registry), `prompts-and-templates.md` (FastMCP GitHub, docs, MCP spec prompts), `typescript-mcp-server.md` (TypeScript SDK, MCP spec, FastMCP docs, Zod). Auth API citations from 2026-02-21 fact-check preserved inline in development-guidelines.md Sources section.
**Description**: The plugin contains over 1200 lines of FastMCP 3.x API documentation derived from a release candidate version, with no source citations. Per CLAUDE.md citation requirements, all factual claims must have cited sources. Auth API claims verified but citations still needed for all documentation.
**Citations**:
- FastMCP auth: `fastmcp/server/auth/__init__.py` lines 8-13, `authorization.py` lines 48, 78, 106 (accessed 2026-02-21)
- FastMCP middleware: `fastmcp/server/middleware/authorization.py` line 51 (accessed 2026-02-21)
**Files**:
- `plugins/fastmcp-creator/skills/fastmcp-creator/references/` (multiple files)

### Consolidate validate_frontmatter.py into plugin_validator.py

**Source**: Frontmatter validation bug-fix session 2026-02-20
**Added**: 2026-02-20
**Completed**: 2026-02-20
**Status**: DONE — commit bfb03f1 on branch claude/fix-frontmatter-validation-u2iV2
**Priority**: P1
**Description**: `plugin_validator.py` (4190 lines) copy-pastes frontmatter
validation logic from `validate_frontmatter.py` (1341 lines) instead of
importing it. Comments acknowledge this: "PYDANTIC FRONTMATTER MODELS (from
validate_frontmatter.py)" and "Complexity preserved from validate_frontmatter.py
for behavioral parity." This creates two maintenance surfaces: any change must
be applied to both scripts (as seen when reversing the name-field bug workaround).

**Required work:**
1. Audit both scripts for all divergences (validation checks, Pydantic models,
   helper functions, CLI flags, scan patterns).
2. Extract shared code (Pydantic models, `_fix_skill_name*`, `extract_frontmatter`,
   `detect_file_type`, `validate_and_normalize`) into a shared module
   `plugins/plugin-creator/scripts/frontmatter_core.py`.
3. Refactor `validate_frontmatter.py` and `plugin_validator.py` to import from
   `frontmatter_core.py`.
4. Port any validation steps present in `validate_frontmatter.py` but missing
   from `plugin_validator.py` (e.g. skill-directory-name check, name-mismatch
   warning added in 2026-02-20).
5. Update `frontmatter_utils.py` if overlap exists.
6. Update all documentation: CLAUDE.md, reference files, script docstrings.
7. Update tests to import from the correct locations.
8. Run full test suite; verify pre-commit hooks still pass.

**Suggested location**: `plugins/plugin-creator/scripts/`

### Validate and verify orchestrator-discipline plugin hooks and processes

**Source**: Plugin creation session 2026-02-19
**Added**: 2026-02-19
**Completed**: 2026-02-21
**Status**: DONE — commits 49e0ae5, ea3e737 on branch claude/bulk-backlog-grooming-LIQDi
**Plan**: plan/tasks-4-validate-orchestrator-discipline.md
**Description**: T1-T4 complete. Plugin passes `claude plugin validate`. Hook directory detection added. `user-invocable: true` added to SKILL.md. All 5 hook behavior tests pass.
**Suggested location**: `plugins/orchestrator-discipline/`


### SAM: Error Recovery / Rollback Procedures

**Source**: Gap analysis of SAM framework
**Added**: 2026-02-01
**Plan**: plan/tasks-5-sam-error-recovery.md
**Description**: Define explicit procedure when a task fails irrecoverably. How to undo artifact changes? How to restore artifact plane to consistent state after failure?
**Research first**: How do GSD, BMAD-METHOD, AutoGPT, and traditional CI/CD handle rollback? What patterns exist for transactional artifact updates?
**Suggested location**: [`stateless-software-engineering-framework.md`](https://github.com/bitflight-devops/stateless-agent-methodology/blob/main/stateless-software-engineering-framework.md) (new Appendix or Part 6 addition)

### SAM: Human Escalation Criteria

**Source**: Gap analysis of SAM framework
**Added**: 2026-02-01
**Description**: Define explicit triggers for escalating to human at each stage (not just Discovery). When should agents block and ask vs attempt repair vs fail?
**Research first**: How do GSD deviation rules work? How does BMAD-METHOD handle human checkpoints? What escalation patterns exist in agent frameworks?
**Suggested location**: [`stateless-software-engineering-framework.md`](https://github.com/bitflight-devops/stateless-agent-methodology/blob/main/stateless-software-engineering-framework.md) (each Agent Specification section)

### SAM: Timeout/Stall Detection

**Source**: Gap analysis of SAM framework
**Added**: 2026-02-01
**Description**: Define mechanism to detect when an agent is stuck or has stalled. Include timeout thresholds per stage, health check patterns, and recovery actions.
**Research first**: How do orchestration frameworks (Temporal, Prefect, Airflow) handle task timeouts? What heartbeat patterns exist? How does Gas Town handle session recycling?
**Suggested location**: [`stateless-software-engineering-framework.md`](https://github.com/bitflight-devops/stateless-agent-methodology/blob/main/stateless-software-engineering-framework.md) (Orchestrator section 3.8)

### SAM: Artifact Schema Validation

**Source**: Gap analysis of SAM framework
**Added**: 2026-02-01
**Description**: Define formal validation rules or JSON schemas for artifact formats. Currently only templates provided. Enable automated validation at stage boundaries.
**Research first**: How do GSD artifacts (STATE.md, ROADMAP.md) enforce structure? What validation approaches exist in BMAD-METHOD? JSON Schema vs YAML validation vs custom parsers?
**Suggested location**: [`sam-artifact-schemas/`](https://github.com/bitflight-devops/stateless-agent-methodology) (new directory with schema files)

### SAM: Scope Creep Detection

**Source**: Gap analysis of SAM framework
**Added**: 2026-02-01
**Description**: Define mechanism to detect when execution diverges from plan. How does Forensic Review detect that the execution agent solved a different problem than planned?
**Research first**: How does GSD plan-checker detect deviation? What diff/comparison techniques exist? How do code review tools detect scope creep in PRs?
**Suggested location**: [`stateless-software-engineering-framework.md`](https://github.com/bitflight-devops/stateless-agent-methodology/blob/main/stateless-software-engineering-framework.md) (section 3.6 Forensic Review)


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
**Target files**: [`stateless-agent-methodology.md`](https://github.com/bitflight-devops/stateless-agent-methodology/blob/main/stateless-agent-methodology.md), [`stateless-software-engineering-framework.md`](https://github.com/bitflight-devops/stateless-agent-methodology/blob/main/stateless-software-engineering-framework.md)
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

### ~~Add ty support alongside mypy in distributed plugins~~ COMPLETED

**Source**: Plugin code review session 2026-02-21
**Added**: 2026-02-21
**Completed**: 2026-02-21
**Description**: Completed as full mypy-to-ty migration across all active documentation, skills, agents, and reference files. All stale mypy references updated to ty (Astral ecosystem). Third-party reference docs (mypy-docs/, rules/mypy/) retained as-is. Historical plan/ files left unchanged.

### conventional-commits: Fix CHANGELOG references to nonexistent files

**Source**: Plugin code review session 2026-02-21
**Added**: 2026-02-21
**Description**: CHANGELOG references files that do not exist in the repository. Additionally, related skills referenced in the plugin do not exist. All dead references need to be either created or removed.
**Files**: `plugins/conventional-commits/` (CHANGELOG and skill cross-references)

### dasel: Reconcile 265 `-f` flag occurrences with reference documentation

**Source**: Plugin code review session 2026-02-21
**Added**: 2026-02-21
**Description**: Reference documentation states dasel uses a specific flag pattern, but 265 occurrences of the `-f` flag across the skill contradict this documentation. The hook file exists on disk but is not registered in the plugin manifest (`plugin.json`). Run `auto_sync_manifests.py --reconcile` to fix manifest drift, and audit `-f` flag usage against official dasel documentation.
**Files**:
- `plugins/dasel/` (skill files with `-f` flag usage)
- `plugins/dasel/.claude-plugin/plugin.json` (missing hook registration)

### ~~litellm: Remove private API documentation and update verification date~~ RESOLVED

**Source**: Plugin code review session 2026-02-21
**Added**: 2026-02-21
**Resolved**: 2026-02-21
**Status**: FACT-CHECKED 2026-02-21 — Review claim REFUTED. `litellm._should_retry()` is NOT a private/internal API. It is explicitly documented in official litellm Exception Mapping docs (<https://docs.litellm.ai/docs/exception_mapping>) with code examples showing `litellm._should_retry(e.status_code)`. Despite the underscore prefix, this is an intentionally public utility function exported in `__init__.py`. The function exists at `utils.py:6506`.
**Description**: Original review flagged `_should_retry()` as private API. Fact-checking verified it is documented public API in official docs. Stale verification date may still need updating.
**Citations**:
- litellm docs: <https://docs.litellm.ai/docs/exception_mapping> (accessed 2026-02-21)
- litellm source: `litellm/utils.py` line 6506 (accessed 2026-02-21)
- litellm source: `litellm/__init__.py` type stub exports `_should_retry` (accessed 2026-02-21)
**Files**:
- `plugins/litellm/skills/litellm/SKILL.md` (line 257)

### verification-gate: Remove unsubstantiated 95% confidence claim

**Source**: Plugin code review session 2026-02-21
**Added**: 2026-02-21
**Description**: The skill contains an unsubstantiated claim about "95% confidence" with no supporting data, citation, or methodology. Per CLAUDE.md verification protocol, claims must be cited. Either add supporting evidence or remove the specific percentage. Also fix missing code fence language specifiers and stale dates.
**Files**: `plugins/verification-gate/` (SKILL.md and reference files)

### development-harness: Remove hardcoded machine path and fix version drift

**Source**: Plugin code review session 2026-02-21
**Added**: 2026-02-21
**Description**: Contains a hardcoded machine-specific path (likely `/home/user/...` or similar) that won't work on other systems. Version references have drifted from actual tool versions. Role table has inconsistencies between documented and actual roles.
**Files**: `plugins/development-harness/` (SKILL.md and reference files)

### llamafile: Fix HuggingFace model URLs (wrong org name + fabricated repos)

**Source**: Plugin code review session 2026-02-21
**Added**: 2026-02-21
**Priority**: P2
**Status**: FACT-CHECKED 2026-02-21 — SourceForge URL VERIFIED (legitimate auto-mirror with explicit disclaimer). HuggingFace model URLs REFUTED — wrong org name (`Mozilla` should be `mozilla-ai`) and model repo names (`gemma-3-3b-it-gguf`, `Qwen3-0.6B-gguf`, `Mistral-7B-gguf`, `Llama-3.1-8B-gguf`) all return 404. GitHub URLs using `mozilla-ai/llamafile` are correct (old `Mozilla-Ocho` redirects).
**Description**: SourceForge mirror is confirmed legitimate (auto-mirror with disclaimer). HuggingFace URLs need two fixes: (1) org name `Mozilla` must be changed to `mozilla-ai`, (2) model repo names need verification against actual `mozilla-ai` org repos — current names appear fabricated. The `llava-v1.5-7b-llamafile` URL works only via redirect.
**Citations**:
- SourceForge mirror: <https://sourceforge.net/projects/llamafile.mirror/files/0.9.3/> — disclaimer states "exact mirror of the llamafile project" (accessed 2026-02-21)
- HuggingFace 404s: `Mozilla/Qwen3-0.6B-gguf`, `Mozilla/Mistral-7B-gguf`, `Mozilla/gemma-3-3b-it-gguf`, `Mozilla/Llama-3.1-8B-gguf` (accessed 2026-02-21)
- HuggingFace redirect: `Mozilla/llava-v1.5-7b-llamafile` redirects to `mozilla-ai/llava-v1.5-7b-llamafile` (accessed 2026-02-21)
- GitHub: `mozilla-ai/llamafile` resolves correctly; `Mozilla-Ocho/llamafile` redirects (accessed 2026-02-21)
**Files**:
- `plugins/llamafile/skills/llamafile/SKILL.md` (lines 78, 87-92 model URLs; line 69 SourceForge OK)

### prompt-optimization: Fix unreachable reference files and raw JSX/MDX markup

**Source**: Plugin code review session 2026-02-21
**Added**: 2026-02-21
**Description**: Two reference files are unreachable (not linked from SKILL.md or any other file). Contains raw Anthropic JSX/MDX markup that should be converted to standard markdown for compatibility with Claude Code's markdown rendering.
**Files**: `plugins/prompt-optimization-claude-45/` (reference files)

### brainstorming-skill: Remove orphaned bibliography entry and cross-reference headings

**Source**: Plugin code review session 2026-02-21
**Added**: 2026-02-21
**Description**: Contains an orphaned bibliography entry (referenced nowhere) and orphaned cross-reference section headings that point to removed or renamed content.
**Files**: `plugins/brainstorming-skill/` (SKILL.md and reference files)

### uv: Fix incorrect script paths in README

**Source**: Plugin code review session 2026-02-21
**Added**: 2026-02-21
**Description**: README contains incorrect script paths that don't match the actual file locations. The README is also disproportionately large for what is essentially a thin wrapper plugin.
**Files**: `plugins/uv/` (README.md)

### plugin-creator: Remove dead code and triplicated regex

**Source**: Plugin code review session 2026-02-21
**Added**: 2026-02-21
**Description**: Contains triplicated regex patterns (same regex defined 3 times), a dead `skipped` list that is populated but never read, an unused `sum()` call, and HK005 warning is incorrectly treated as an error in certain code paths. Also has a `noqa BLE001` suppression that should be addressed per CLAUDE.md linting policy.
**Files**: `plugins/plugin-creator/` (scripts and skill files)

### Add PR003/PR004 test coverage to plugin registration validator

**Source**: Code review session 2026-02-21
**Added**: 2026-02-21
**Description**: `PluginRegistrationValidator` defines PR003 (missing metadata fields: repository, homepage, author) and PR004 (repository URL mismatches git remote URL) at lines 276-277 of `plugin_validator.py`, and emits them at lines 2815 and 2834. Tests exist for PR001 (unregistered) and PR002 (missing file), but not PR003/PR004. Add tests to `plugins/plugin-creator/tests/test_plugin_registration_validator.py` covering: (1) PR003 emitted when metadata fields absent; (2) PR004 emitted when repo URL mismatches remote.

---

### kaizen: MCP consolidation analysis

**Source**: Design session 2026-02-20
**Added**: 2026-02-20
**Description**: The plugin currently runs two MCP servers (`kaizen-duckdb` via mcp-server-motherduck, `kaizen-analysis` via server.py) plus a standalone CLI script (`sentiment-score.py`). Investigate: (1) What does each MCP server provide that the other cannot? Can they be merged into a single server? (2) Why is `sentiment-score.py` a standalone script rather than an MCP tool inside `server.py`? What would be gained or lost by moving scoring into the MCP server (always-on scoring, no manual invocation, lock ownership)? (3) Is there a clean boundary between "batch processing" (script) and "query/serve" (MCP) that should be preserved?
**Decision needed**: Consolidate vs. keep separate, with rationale.
**Suggested location**: `plugins/agentskill-kaizen/`

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
**Suggested location**: [`stateless-software-engineering-framework.md`](https://github.com/bitflight-devops/stateless-agent-methodology/blob/main/stateless-software-engineering-framework.md) (new section 2.4 or Appendix)

### SAM: Multi-Model Strategy

**Source**: Gap analysis of SAM framework
**Added**: 2026-02-01
**Description**: Define guidance for using different models for different agent types. E.g., cheaper/faster models for simple verification, stronger models for planning.
**Research first**: How do agent frameworks handle model selection? What cost/quality tradeoffs exist? How does Claude Code's haiku/sonnet/opus selection work?
**Suggested location**: [`stateless-software-engineering-framework.md`](https://github.com/bitflight-devops/stateless-agent-methodology/blob/main/stateless-software-engineering-framework.md) (Implementation Roadmap or new Appendix)

### SAM: Audit Trail / Observability

**Source**: Gap analysis of SAM framework
**Added**: 2026-02-01
**Description**: Beyond artifacts, define logging/metrics/tracing guidance. How to diagnose pipeline issues? What telemetry to capture?
**Research first**: How do GSD and BMAD-METHOD handle logging? What observability patterns exist in agent frameworks? OpenTelemetry for LLM workflows?
**Suggested location**: [`stateless-software-engineering-framework.md`](https://github.com/bitflight-devops/stateless-agent-methodology/blob/main/stateless-software-engineering-framework.md) (new Appendix I)

### SAM: Partial Success Handling

**Source**: Gap analysis of SAM framework
**Added**: 2026-02-01
**Description**: Define how to represent and handle partial task success. Task completes some DoD items but not all. How is this state represented in artifacts?
**Research first**: How do GSD checkpoints represent partial progress? How do CI/CD systems handle partial test passes? What state machine patterns exist?
**Suggested location**: [`stateless-software-engineering-framework.md`](https://github.com/bitflight-devops/stateless-agent-methodology/blob/main/stateless-software-engineering-framework.md) (section 3.5 Execution Agent output)

### SAM: Context Size Management

**Source**: Gap analysis of SAM framework
**Added**: 2026-02-01
**Description**: Define explicit guidance for measuring and managing context size per agent. What's the target token budget? How to detect context pressure?
**Research first**: How do agent frameworks measure context usage? What token counting approaches exist? How does Claude Code handle context limits internally?
**Suggested location**: [`stateless-software-engineering-framework.md`](https://github.com/bitflight-devops/stateless-agent-methodology/blob/main/stateless-software-engineering-framework.md) (section 2.1 or Appendix C)

### SAM: Conflicting Review Findings

**Source**: Gap analysis of SAM framework
**Added**: 2026-02-01
**Description**: Define protocol when forensic review and self-verification disagree. Which takes precedence? How to adjudicate conflicts?
**Research first**: How do code review systems handle conflicting reviewers? What adjudication patterns exist in multi-agent systems? How does GSD handle verification disagreements?
**Suggested location**: [`stateless-software-engineering-framework.md`](https://github.com/bitflight-devops/stateless-agent-methodology/blob/main/stateless-software-engineering-framework.md) (section 3.6 Forensic Review)

### Multi-session build state lost during context compaction

**Source**: agentskill-kaizen plugin build (2026-02-18), 3 sessions with 1 compaction
**Added**: 2026-02-18
**Description**: During the agentskill-kaizen build (8-phase `/plugin-dev:create-plugin` workflow), context compaction mid-build converted structured task state into a narrative summary. The resuming session had to reconstruct "what's done vs pending" from prose rather than a checklist. Background agent results that were already consumed and applied reappeared as late notifications after compaction, requiring manual deduplication ("did I already handle this?"). No persistent artifact tracked phase completion, commit SHAs per phase, deferred items, or agent result consumption status.
**Observed symptoms**:
- Phase completion status existed only in ephemeral context — lost on compaction
- Background agent notifications arrived after their findings were already applied (3 duplicate notifications)
- Plan committed early (`87a0b93`) diverged from actual implementation but was never updated
- No mechanism to mark agent results as "consumed" — each notification required re-evaluation

### `/plugin-dev:create-plugin` workflow lacks intra-phase parallelism tracking

**Source**: agentskill-kaizen plugin build (2026-02-18)
**Added**: 2026-02-18
**Description**: The 8-phase create-plugin workflow treats each phase as a serial step, but Phase 5 (Implementation) actually consisted of 6 parallel sub-tasks and Phase 6 (Validation) spawned 3 parallel review agents. The workflow provides no structure for tracking parallel work within a phase — no task dependencies, no completion gates, no way to know which sub-tasks are done after compaction. Batching validation fixes by file rather than by finding would also have been more efficient (SKILL.md was edited 3 separate times when one pass would have sufficed).

### Background agent result deduplication after compaction

**Source**: agentskill-kaizen plugin build sessions 2-3 (2026-02-18)
**Added**: 2026-02-18
**Description**: Background agents launched in session N may complete after context compaction or session restart. The system delivers their results as `<task-notification>` messages, but there is no mechanism to mark results as already consumed. During the kaizen build, 3 review agents (plugin-validator, 2x skill-reviewer) completed during Phase 6 and their findings were applied in commit `0d61480`. After compaction, all 3 re-delivered their notifications in session 3, requiring manual evaluation each time ("was this already handled?"). A persistent state file (e.g., `.planning/kaizen/agent-results.json` tracking task IDs → consumed/pending) would eliminate this waste.

### `/plugin-dev:create-plugin` Phase 6 validation should batch fixes by file, not by finding

**Source**: agentskill-kaizen plugin build Phase 6 (2026-02-18)
**Added**: 2026-02-18
**Description**: Phase 6 collected findings from 3 parallel review agents, then applied fixes one finding at a time. This resulted in SKILL.md being edited 3 separate times (description rewrite, SQL removal, MCP server name fix) when a single pass through the file would have applied all fixes together. The workflow should group all findings by file, then make one editing pass per file. Reduces Edit tool calls and context consumed by repeated reads.

### Plan artifact diverges from implementation without update mechanism

**Source**: agentskill-kaizen plugin build (2026-02-18), plan committed as `87a0b93`
**Added**: 2026-02-18
**Description**: The research/plan from Phase 2-3 was committed early as a markdown file. During implementation (Phase 5), decisions changed — the MCP server grew from planned scope, analysis dimensions were rebalanced, hook patterns shifted. The plan was never updated to reflect actual implementation. After compaction, the stale plan became a potential source of confusion. Two options to address: (1) update the plan artifact after each phase, or (2) treat the plan as disposable and track only the living state (what's done, what's pending, what deferred).

### Evaluate scikit-learn dependency weight for agentskill-kaizen cluster_sessions tool

**Source**: agentskill-kaizen MCP server review (2026-02-18)
**Added**: 2026-02-18
**Description**: The `cluster_sessions` tool in `plugins/agentskill-kaizen/mcp/server.py` uses `scikit-learn` (`KMeans`, `CountVectorizer`) for session clustering. scikit-learn pulls in ~40MB of transitive dependencies (numpy, scipy, joblib, threadpoolctl). The typical use case is clustering dozens of sessions, not thousands. The Phase 1-2 research produced no durable artifact evaluating library choices — `scikit-learn` was assumed from the backlog without validation. Investigate whether a lighter alternative (e.g., `pyclustering`, stdlib-based implementation, or just `numpy` directly) would suffice for this scale, and whether the KMeans-on-bag-of-words approach is even appropriate for tool-call sequence similarity.
**Research first**: What clustering approaches work for short categorical sequences? Is cosine similarity on bag-of-words tool vectors meaningful for workflow comparison? What lightweight Python clustering libraries exist that don't pull in scipy?

---

### github_project_setup.py: add `milestone close` command

**Source**: PR #149 follow-up — start-milestone automation (2026-02-22)
**Added**: 2026-02-22
**Completed**: 2026-02-22
**Status**: DONE — Added `milestone close` subcommand to `github_project_setup.py`: validates milestone is open, lists open issues, transitions `status:in-progress` → `status:done`, closes milestone, prints summary. Added `status:done` to LABELS taxonomy. Added `_transition_to_done()` helper. Updated `complete-milestone/SKILL.md` Step 4 to reference the script command (consistent with `start-milestone`).
**Description**: `github_project_setup.py` now has `milestone start` which bulk-transitions `status:needs-grooming` → `status:in-progress`. Add a symmetric `milestone close` command for the `complete-milestone` skill. The command should: (1) validate the milestone is open, (2) list all still-open issues (warn if any remain), (3) transition open issues from `status:in-progress` → `status:done` or close them, (4) close the milestone itself via `milestone.edit(state="closed")`, (5) print a completion summary. Update `complete-milestone/SKILL.md` to reference the script (consistent with `start-milestone`).
**Suggested location**: `.claude/skills/gh/scripts/github_project_setup.py` — add `milestone close` subcommand under `milestone_app`; update `.claude/skills/complete-milestone/SKILL.md`

---

### work-backlog-item: accept `#N` in `close` and `resolve` routing

**Source**: PR #149 follow-up — start-milestone automation (2026-02-22)
**Added**: 2026-02-22
**Completed**: 2026-02-22
**Status**: DONE — Updated `work-backlog-item/SKILL.md`: Arguments table now shows `close #N` / `resolve #N` forms; Routing table updated with `(title or #N)` annotation; Step 9a now handles `#N` argument by fetching title from GitHub Issue via `gh issue view`; Step 9 Extension updated to use issue number from `$1` when invoked as `close #N`; Error Handling and Example Sessions updated.
**Description**: `work-backlog-item close` currently expects a title substring. Since GitHub Issues are now the canonical source of truth, `close #N` and `resolve #N` should also be valid — fetch the title from the issue and proceed. This eliminates the need to remember exact title substrings for the common case of closing work started with `/work-backlog-item #N`. Routing table in `SKILL.md` needs a new `#N` row for the close/resolve paths.
**Suggested location**: `.claude/skills/work-backlog-item/SKILL.md` — Routing table (Step 9 path column), plus Step 9 Extension GitHub issue close command

---

### github_project_setup.py: add GitHub Projects V2 status field updates

**Source**: PR #149 follow-up — start-milestone automation (2026-02-22)
**Added**: 2026-02-22
**Description**: The `milestone start` and `milestone close` commands transition `status:*` labels on issues but do not update the GitHub Projects V2 `Status` field (the kanban column). These are separate systems — labels are on the issue, the V2 status is a project-level field. Add `project update-status` subcommand to `github_project_setup.py` that uses the GraphQL API to set the `Status` field on project items. This unlocks the full kanban board view in GitHub Projects.
**Research first**: Projects V2 GraphQL API — `updateProjectV2ItemFieldValue` mutation. Field ID discovery via `projectV2Fields` query. Reference: `.claude/skills/gh/references/projects-v2.md`.
**Suggested location**: `.claude/skills/gh/scripts/github_project_setup.py` — add `project_app` Typer sub-app with `project update-status --project-number N --issue-number N --status "In Progress"` command

---

### Configurable Token Thresholds for plugin_validator.py

**Source**: Session 2026-02-18, token threshold analysis
**Added**: 2026-02-18
**Description**: `plugin_validator.py` token complexity thresholds (4400 warn, 8800 error) are hardcoded constants. Allow per-project or per-plugin configuration via a config file (e.g., `.claude-plugin/validator.json` or a `[tool.plugin-validator]` section in `pyproject.toml`) so different projects can set their own limits.
**Suggested location**: `plugins/plugin-creator/scripts/plugin_validator.py`

---

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

### Resolve 48 pre-existing ty (Astral type checker) diagnostics

**Source**: `uv run ty check .` run during session 2026-02-11
**Completed**: 2026-02-19
**Description**: ty v0.0.16 reported 48 diagnostics. By v0.0.17 with existing `extra-paths` config, count dropped to 34 warnings (0 errors), all in `plugins/agentskill-kaizen/tests/`. Root causes: FastMCP `@mcp.tool()` decorator returns `FunctionTool` (ty sees as non-callable), and `ModuleType` dynamic attribute stubs. Suppressed `call-non-callable` and `unresolved-attribute` to `"ignore"` in test-file override scope. Added `typecheck-ty` CI job in `allowed-failures` for conservative rollout.
**Location**: `pyproject.toml` (`[tool.ty.overrides]`), `.github/workflows/code-quality.yml`

### Resolve plugin-validator pre-existing errors to make CI gate blocking

**Source**: CI pipeline run 22018867027 on claude/fix-ci-pipeline-RJ0Tw (2026-02-14)
**Completed**: 2026-02-18
**Description**: Original 19 errors and 63 warnings (from 2026-02-14) were resolved across multiple prior sessions: LK001 broken links fixed, SK007 orchestrating-swarms split into 4 sub-skills, SK005 trigger phrase issues resolved. Final session fixed 4 remaining LK002 warnings in agentskills/SKILL.md (added `./` prefix) and promoted `validate-plugins` from `allowed-failures` to blocking in the quality gate.
**Location**: `plugins/plugin-creator/skills/agentskills/SKILL.md`, `.github/workflows/code-quality.yml`
**Plan**: plan/tasks-3-validator-ci-gate.md

### Resolve manifest drift to make CI sync check blocking

**Source**: CI pipeline run 22018867027 on claude/fix-ci-pipeline-RJ0Tw (2026-02-14)
**Completed**: 2026-02-18
**Description**: Ran `auto_sync_manifests.py --reconcile` to sync all 22 plugin.json manifests with on-disk skills, commands, and agents (153 drift entries fixed). Split the `validate-plugins` CI job into two: `validate-plugins` (plugin-validator, stays in `allowed-failures`) and `manifest-sync` (now blocking in quality gate). Version bumps applied to all affected plugins.
**Location**: `plugins/*/. claude-plugin/plugin.json`, `.github/workflows/code-quality.yml`

### Replace requests with httpx in all scripts

**Source**: CI/pre-commit inconsistency discovery (2026-02-05)
**Completed**: 2026-02-06
**Description**: Migrated `validate_glfm.py` from `requests` to `httpx`. Added TID251 ruff ban rule for `requests` imports. Removed `types-requests` from dev dependencies. `sync_gitlab_docs.py` was already using httpx.
**Location**: `plugins/gitlab-skill/skills/gitlab-skill/scripts/validate_glfm.py`, `pyproject.toml`

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
**Location**: [`stateless-software-engineering-framework.md`](https://github.com/bitflight-devops/stateless-agent-methodology/blob/main/stateless-software-engineering-framework.md) (section 2.1.2)

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

**Skill for research phase**: `/research-and-compare <methodology-name-or-url>` (moved to [stateless-agent-methodology](https://github.com/bitflight-devops/stateless-agent-methodology) repo)

- Produces structured comparison documents following SAM comparison template
- Includes overlap/divergence analysis, weakness discovery, implementation pairing
- Outputs to [`stateless-agent-methodology/.meta/v1_comparisons/`](https://github.com/bitflight-devops/stateless-agent-methodology/tree/main/.meta/v1_comparisons)

**Existing SAM comparisons** (start here before running new research):

- [SAM v1 comparisons](https://github.com/bitflight-devops/stateless-agent-methodology/tree/main/.meta/v1_comparisons)
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
2. If more research needed: `/research-and-compare <framework>` (in stateless-agent-methodology repo) for specific topics
3. Synthesize findings into SAM framework update
4. Mark backlog item complete

### P1: plugin-validator pre-commit output is too noisy

**Source**: `prek --all-files` run during session 2026-02-15
**Added**: 2026-02-15
**Description**: Running `prek --all-files` produces pages of plugin-validator output — every file in the repo gets validated, and pre-existing warnings/errors (broken links, SK006 complexity, missing trigger phrases) for files unrelated to the current change drown out actionable information. The validator needs a review of its pre-commit integration to reduce noise.
**Observed symptoms**:
- Full `prek --all-files` run produces hundreds of lines of validator output
- Pre-existing LK001 (broken links), SK006 (complexity warnings), SK005 (missing triggers) reported on every run
- FM003 (no frontmatter) fires on template files (`commands/development/templates/*.md`) that are not components — they're templates meant to be copied and filled in, not commands/skills/agents
- Difficult to find actual issues in changed files among the noise
**Related**: "plugin-validator UX and coverage gaps" sub-issues below, "Resolve plugin-validator pre-existing errors to make CI gate blocking" backlog item
**File**: `plugins/plugin-creator/scripts/plugin_validator.py`, `.pre-commit-config.yaml`

### P1: plugin-validator UX and coverage gaps

**Source**: Experimental validation of plugin_validator.py against all component types (2026-02-13)
**Added**: 2026-02-13
**Last groomed**: 2026-02-13
**Plan**: plan/tasks-2-validator-ux-coverage.md
**File**: `./plugins/plugin-creator/scripts/plugin_validator.py` (2934+ lines)
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

Use `/python3-development:add-new-feature` to extend `plugin_validator.py`. Delegate implementation to `@python-cli-architect`. Each sub-issue is a separate task.
