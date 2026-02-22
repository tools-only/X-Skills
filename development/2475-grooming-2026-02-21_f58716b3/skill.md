# Backlog Grooming Report

**Date**: 2026-02-21
**Items groomed**: 39 active (1 DONE skipped)
**Arguments**: all
**Groomer agents launched**: 10 (2 waves of 5) + 1 single-item
**Session branch**: claude/bulk-backlog-grooming-LIQDi

---

## Summary

| Item | Priority | RT-ICA | Research Found | Blockers | Effort |
|------|----------|--------|----------------|----------|--------|
| hallucination-detector: Fix dead backtick evidence marker and multi-occurrence "because" bug | P1 | APPROVED | Hook source, plugin files | Missing tests | Small |
| plugin-validator UX and coverage gaps | P1 | APPROVED | Skills, agents, test files | 0 | Medium |
| Validate orchestrator-discipline plugin | P1 | APPROVED | Plan file, hook files | 0 | Medium |
| SAM Extension — Integrate ARL General Theory | P1 | APPROVED | ARL references, 7 principles | 0 | Medium (cross-repo) |
| ARL Skill Development | P1 | APPROVED | ARL scaffold exists | Design question | Medium |
| Meta-Process Capture — Expert Panel Dataset Builder | P1 | APPROVED | Command exists | Design decision | Medium |
| SAM: Error Recovery / Rollback Procedures | P1 | APPROVED | Comparison docs | Research-first + cross-repo | High |
| SAM: Human Escalation Criteria | P1 | APPROVED | Comparison docs | Research-first + cross-repo | High |
| SAM: Timeout/Stall Detection | P1 | APPROVED | Comparison docs | Research-first + cross-repo | High |
| SAM: Artifact Schema Validation | P1 | APPROVED | Comparison docs | Research-first + cross-repo | High |
| SAM: Scope Creep Detection | P1 | APPROVED | Comparison docs | Research-first + cross-repo | High |
| Extract claude-plugin-lint to standalone PyPI package | P1 | APPROVED | Source file | Missing packaging setup | Very High |
| plugin-validator pre-commit noise | P1 | APPROVED | CLI flags exist (unused) | Design option needed | Small-Medium |
| Add PR003/PR004 test coverage | P2 | APPROVED | Test patterns (PR001/PR002) | 0 | Small |
| SAM: Parser regex false positive | P2 | APPROVED | File path, line number | No unit tests | Small |
| kaizen: MCP consolidation analysis | P2 | APPROVED | Plugin location known | None | Small-Medium |
| SAM: Replace validate-task-file.sh | P2 | APPROVED | task_format.py found | Unclear invocation | Small |
| Multi-session build state lost | P2 | APPROVED | STATE.md pattern found | Design work | Medium |
| /plugin-dev intra-phase parallelism | P2 | APPROVED | Workflow location found | Design needed | Medium |
| Background agent result deduplication | P2 | APPROVED | Solution proposed | None | Small |
| /plugin-dev Phase 6 batch fixes | P2 | APPROVED | Workflow found | None | Small |
| Plan artifact diverges from implementation | P2 | APPROVED | Two options stated | Option selection needed | Small |
| Evaluate scikit-learn dependency weight | P2 | APPROVED | Plugin location known | None | Small |
| SAM: Parallel Execution Details | P2 | APPROVED | Comparison docs | Research-first + cross-repo | High |
| SAM: Multi-Model Strategy | P2 | APPROVED | Claude docs | Research-first + cross-repo | Medium |
| SAM: Audit Trail / Observability | P2 | APPROVED | Comparison docs | Research-first + cross-repo | High |
| SAM: Partial Success Handling | P2 | APPROVED | Comparison docs | Research-first + cross-repo | High |
| SAM: Context Size Management | P2 | APPROVED | Claude docs | Research-first + cross-repo | Medium |
| SAM: Conflicting Review Findings | P2 | APPROVED | Comparison docs | Research-first + cross-repo | Medium |
| Configurable Token Thresholds | Idea | APPROVED | plugin_validator.py | None | Small — promote to P2 |
| SAM: Cost/Token Management | Idea | APPROVED | — | Research-first | Explore |
| SAM: Team Coordination Protocols | Idea | APPROVED | — | Research-first | Explore |
| SAM: External System Integration | Idea | APPROVED | — | Research-first | Explore |
| SAM: Migration Strategy Guide | Idea | APPROVED | — | Research-first | Explore |
| SAM: Training/Onboarding Materials | Idea | APPROVED | — | Research-first | Explore |
| SAM: Non-Code Workflow Guidance | Idea | APPROVED | — | Research-first | Explore |
| Carbonyl Browser Integration | Idea | BLOCKED | — | Needs TTY + network | Hold |
| Validate is-fast | Idea | BLOCKED | — | DNS blocked | Hold |
| Validate agent-browser | Idea | BLOCKED | — | Playwright download blocked | Hold |

---

## Individual Manifests

### hallucination-detector: Fix dead backtick evidence marker and multi-occurrence "because" bug (P1)

**RT-ICA**: APPROVED — both bugs VERIFIED from source, 1 assumption (test file location)

**Fact-check**: 2/2 VERIFIED from direct source read (primary source: `plugins/hallucination-detector/scripts/hallucination-audit-stop.js`)

**Bug 1 (lines 111-114 vs 127-141)**: `stripLowSignalRegions()` removes inline code spans before `hasEvidenceNearby()` is called. `EVIDENCE_MARKERS[0]` (`/\`[^\`]+\`/`) can never match in the already-stripped haystack. Fix: pass original unstripped text to`hasEvidenceNearby()`.

**Bug 2 (lines 223-229, 306-332)**: All phrase loops use `lower.indexOf(phrase)` which returns only the first occurrence index. Subsequent occurrences of the same phrase are never evaluated. Fix: replace `indexOf` with `while`/`matchAll` loop.

**Skills**: none specifically — standard JS fix
**Agents**: @javascript-pro (matchAll pattern), @code-review
**Prior work**: `plugins/hallucination-detector/README.md`, `plugins/hallucination-detector/.claude-plugin/plugin.json`
**Missing**: No test file for `hallucination-audit-stop.js` — create test fixtures in `plugins/hallucination-detector/tests/`
**Effort**: Small (2 focused logic changes in one file)

**Suggested first steps**:
1. Create test fixtures with backtick-evidence + multi-occurrence "because" scenarios
2. Fix Bug 1: thread original `text` into `hasEvidenceNearby()` calls alongside stripped `haystack`
3. Fix Bug 2: replace `indexOf()` with `matchAll()` loop in speculationPhrases and causalityPhrases loops
4. Run `uv run prek run --files plugins/hallucination-detector/scripts/hallucination-audit-stop.js`
5. Validate plugin: `uv run plugins/plugin-creator/scripts/plugin_validator.py plugins/hallucination-detector/`

---

### plugin-validator UX and coverage gaps (P1)

**RT-ICA**: APPROVED — 0 blockers

**Key discoveries (groomer agent ae49af7)**:
- `FileType` enum already has `HOOK_CONFIG` variant (lines 318-329) — not missing as backlog implied
- `test_hook_validator.py` already exists — sub-issue 3 (hook support) has test infrastructure ready
- `DescriptionValidator` already has `file_type` awareness at line 1981 — sub-issue 2 smaller than expected
- 19 test files available (backlog stated 12)
- QA report at `plugins/plugin-creator/planning/plugin-validator-qa-report.md` maps all 4 sub-issues

**Recommended order**: Sub-issue 4 (dead code removal, lines 904-911) → Sub-issue 1 (UX, lines 2928-2931) → Sub-issue 3 (hook support) → Sub-issue 2 (command SK005)

**Skills**: plugin-creator:claude-hooks-reference-2026, python3-development
**Agents**: @refactor-executor, @plugin-assessor
**Prior work**: plan/tasks-2-validator-ux-coverage.md, planning/plugin-validator-qa-report.md

---

### Validate and verify orchestrator-discipline plugin (P1)

**RT-ICA**: APPROVED — 0 external blockers; internal blockers are planned tasks

**Key discoveries (groomer agent a37eae0)**:
- Plugin currently FAILS `claude plugin validate` due to invalid `rules` field in plugin.json
- `SKILL.md` missing `user-invocable: true` frontmatter — skill activation broken
- Hook gap: directory-targeted Grep calls (e.g., `path: "src/"`) are NOT detected
- Plan at `plan/tasks-4-validate-orchestrator-discipline.md` is comprehensive with 4 tasks (T1-T4), architectural decisions pre-made, 5 verification test commands specified

**Task execution order**: T1 (fix plugin.json invalid `rules` field) + T2 (add directory detection to hook) + T3 (add `user-invocable: true`) → T4 (convergence validation)

**Skills**: plugin-creator:claude-plugins-reference-2026
**Agents**: @general-purpose (per plan task assignments)
**Prior work**: plan/tasks-4-validate-orchestrator-discipline.md

---

### Add PR003/PR004 test coverage (P2)

**RT-ICA**: APPROVED — 0 blockers

**Key discoveries (groomer agent abe9d6c)**:
- PR003 = info severity (not warning/error) — tests must check `result.info` list
- PR004 = warning severity — tests must check `result.warnings` list
- Helper functions `_make_plugin()` (lines 32-64) provide complete test setup pattern
- Existing PR001/PR002 test classes (lines 267-543) serve as templates
- Git metadata extraction may require `monkeypatch` for PR004 tests

**Target file**: `plugins/plugin-creator/tests/test_plugin_registration_validator.py`
**New classes**: `TestMissingMetadata` (PR003) and `TestRepositoryMismatch` (PR004)

---

### SAM Extension — Integrate ARL General Theory (P1)

**RT-ICA**: APPROVED — 0 blockers

**Key discoveries (groomer agent af51cae)**:
- **Path correction**: ARL references are at `plugins/plugin-creator/skills/arl/references/` NOT `skills/assessor/references/ARL/` (backlog path is wrong)
- All 7 principles fully documented in `synthesis-general-theory.md` (173 lines) with evidence from 6 frameworks
- `synthesis-arl-applicable.md` provides mapping template (R1-R10 → framework patterns)
- Target files are in external repo `bitflight-devops/stateless-agent-methodology` (confirmed)
- Blocks 5 related P1 SAM gap items

**Requires**: Cross-repo work in `bitflight-devops/stateless-agent-methodology`

---

### ARL Skill Development (P1)

**RT-ICA**: APPROVED — 0 blockers; design question before starting

**Key discoveries (groomer agent a37be14)**:
- ARL SKILL.md scaffold already exists at `plugins/plugin-creator/skills/arl/SKILL.md`
- R1-R10 gates fully specified in `synthesis-arl-applicable.md` (305 lines)
- Design question: **Reference skill** (documents gates) vs. **Orchestration skill** (executes ARL loop)?
- Backlog scope says "logical process design" → reference skill, but assessor SKILL.md shows orchestration pattern
- Soft dependency on SAM Extension (not blocking)

---

### Meta-Process Capture — Expert Panel Dataset Builder (P1)

**RT-ICA**: APPROVED — design decision needed before starting

**Key discoveries (groomer agent a7d0457)**:
- `.claude/commands/arl-expert-panel.md` already exists (100 lines) — process is operational
- `ARL-agent-instructions.md` at `plugins/plugin-creator/skills/arl/references/` (440 lines) — complete process specification
- **Design decision**: Skill format (`.claude/skills/expert-panel-methodology/SKILL.md`) vs. reference document
- Blocks SAM Extension and ARL Skill Development items

---

### SAM: Parser regex false positive (P2)

**RT-ICA**: APPROVED — soft blockers

**Key discoveries (groomer agent aea3b2f)**:
- File: `plugins/python3-development/skills/implementation-manager/scripts/implementation_manager.py` line 645
- Regex: `^#{2,3}\s+Task:?\s+([A-Za-z0-9.]+)[:\s-]+(.+)$`
- Fix: All valid task IDs contain at least one digit (1.1, T1, T15); "Summary" is all-alpha → constrain task ID to require digit
- No existing unit tests for `parse_task_content()` — tests must be added with fix
- `task_status_hook.py` may also be affected

---

### Multi-session build state lost (P2)

**RT-ICA**: APPROVED — 0 blockers

**Key discoveries (groomer agent a8c5cfe)**:
- `STATE.md` artifact pattern already defined in plugin-creator SKILL.md (lines 41-83) — foundation exists
- Solution: Extend STATE.md with phase completion checklist + commit SHA per phase
- Add `agent-results.json` tracking task ID → consumed status
- Addresses all 4 related backlog items (intra-phase parallelism, result deduplication, batch fixes, plan divergence)
- No code changes needed — documentation/artifact format update only

---

### SAM: Replace validate-task-file.sh (P2)

**RT-ICA**: APPROVED — soft blockers

**Key discoveries (groomer agent ad62001)**:
- `validate-task-file.sh` (170 lines) validates old markdown `tasks-refactor-*.md` format
- `task_format.py` exists in `python3-development` plugin (cross-plugin dependency)
- `TASK_FILE_FORMAT.md` at `.claude/docs/TASK_FILE_FORMAT.md` documents target YAML frontmatter format
- **Critical**: Script not found invoked anywhere — may be unused; verify before implementing replacement
- Cross-plugin dependency decision needed for `task_format.py`

---

### plugin-validator pre-commit output is too noisy (P1)

**RT-ICA**: APPROVED — solution approach needs decision

**Key discoveries (groomer agent a5b16d6)**:
- `--filter` and `--filter-type` flags already exist in plugin_validator.py but NOT used in pre-commit hook
- `commands/development/templates/` directory does NOT exist — FM003 complaint may be stale/resolved
- 4 solution options identified:
  - A: Git staging integration (only validate staged files)
  - B: Severity-level filtering via CLI flag
  - C: Project-level config file for exclusions
  - D: `--exclude-codes` CLI flag in hook config
- Lowest-effort quick win: Add `--show-progress` to hook config and use existing `--filter-type` flag

---

## RT-ICA Results

### BLOCKED Items
None — all active items are APPROVED. Three Ideas are blocked by environment constraints:
- Carbonyl Browser Integration: Needs TTY + network
- Validate is-fast: DNS blocked in restricted environment
- Validate agent-browser: Playwright browser download blocked

### APPROVED Items

**P1 (12 active items)**:
- plugin-validator UX gaps: 8 conditions verified
- Validate orchestrator-discipline: 9 conditions verified, 1 missing (test methodology — DERIVABLE)
- SAM Extension — Integrate ARL General Theory: 7 conditions verified
- ARL Skill Development: 7 conditions verified
- Meta-Process Capture — Expert Panel: 7 conditions verified
- SAM: Error Recovery / Rollback Procedures: 4 conditions, research-first
- SAM: Human Escalation Criteria: 4 conditions, research-first
- SAM: Timeout/Stall Detection: 4 conditions, research-first
- SAM: Artifact Schema Validation: 4 conditions, research-first
- SAM: Scope Creep Detection: 4 conditions, research-first
- Extract claude-plugin-lint to PyPI: 5 conditions, missing packaging infrastructure
- plugin-validator pre-commit noise: 8 conditions verified, solution approach pending

**P2 (16 items)**: All APPROVED; 8 with research-first dependencies, 8 immediately actionable

---

## Cross-Item Findings

### Shared Dependencies

1. **plugin_validator.py** — Central dependency for 4 backlog items:
   - plugin-validator UX gaps (direct implementation target)
   - plugin-validator pre-commit noise (configuration target)
   - Add PR003/PR004 test coverage (test target)
   - Configurable Token Thresholds (extension target)
   - **Caution**: Coordinate changes; avoid simultaneous edits

2. **ARL reference files** (`plugins/plugin-creator/skills/arl/references/`):
   - SAM Extension — uses synthesis-general-theory.md
   - ARL Skill Development — uses synthesis-arl-applicable.md
   - Meta-Process Capture — uses ARL-agent-instructions.md, qa-expert-panel.md
   - **Path correction**: Backlog lists `skills/assessor/references/ARL/` — actual path is `skills/arl/references/`

3. **`bitflight-devops/stateless-agent-methodology` repo** — Target for all SAM gap items:
   - All 11 SAM gap items (P1 and P2) require cross-repo access
   - Must clone or checkout repo before any SAM implementation work
   - Confirmed by user: SAM has been moved to its own repo

4. **`/plugin-dev:create-plugin` workflow** — Source of 4 related P2 items:
   - Multi-session build state, intra-phase parallelism, result deduplication, plan divergence
   - Multi-session state item is the root fix; other 3 are corollaries
   - Address together in one session rather than individually

### Suggested Groupings

| Group | Items | Rationale |
|-------|-------|-----------|
| Plugin Validator Fix Bundle | UX gaps + pre-commit noise + PR003/PR004 tests | Same file, coordinate edits |
| ARL Documentation Bundle | SAM Extension + ARL Skill + Expert Panel Meta-Process | Same source files, sequential dependency |
| SAM Gap Bundle | Error Recovery + Human Escalation + Timeout + Artifact Schema + Scope Creep | All cross-repo, research-first, batch together |
| Plugin-Dev Workflow Bundle | Multi-session state + intra-phase parallelism + result dedup + batch fixes + plan divergence | All /plugin-dev:create-plugin workflow, address together |

### Research Gaps

Items requiring external research before implementation (all SAM-related):
- GSD rollback/error handling patterns
- BMAD-METHOD human checkpoint patterns
- Temporal/Prefect/Airflow timeout and heartbeat patterns
- AutoGPT task failure handling
- OpenTelemetry for LLM workflow observability

Recommended methodology: `/research-and-compare <framework>` (in stateless-agent-methodology repo) per backlog Format Guide.

### Pre-Existing Issues Found

1. **Backlog path error**: All ARL-related items reference `plugins/plugin-creator/skills/assessor/references/ARL/` — actual path is `plugins/plugin-creator/skills/arl/references/`. Backlog should be updated.

2. **orchestrator-discipline plugin broken**: Plugin currently fails `claude plugin validate` due to invalid `rules` field. This is a pre-existing functional defect — the plugin is not enforcing delegation discipline as intended.

3. **SKILL.md missing `user-invocable: true`**: orchestrator-discipline skill activation is broken until T3 of plan/tasks-4-validate-orchestrator-discipline.md is executed.

4. **`commands/development/templates/` FM003 complaint**: Directory does not exist — pre-commit noise complaint may be partially stale.

5. **`validate-task-file.sh` not integrated**: Script exists but is not invoked anywhere — unclear if it serves any active purpose.

---

## Top Three Selected Items

### Selection Criteria
1. Priority (P1 > P2)
2. Zero external blockers (fully actionable today)
3. Self-contained in this repository (no cross-repo work)
4. Immediate impact on tooling or security posture
5. Executable by sub-agent with information already gathered

### #1: plugin-validator UX and coverage gaps (P1)

**Rationale**: Highest-impact immediately actionable item. 4 independent sub-issues, all with specific file paths, line numbers, and acceptance criteria. Fixes improve daily developer workflow (false positives suppressed, hooks validated, accurate reporting). No blockers. Plan exists. 19 tests available.

**Dependencies**: None
**Prerequisites**: None
**Complexity**: Medium — Python code changes to plugin_validator.py

**Sub-agent**: `general-purpose`
**Why**: Python implementation task requiring reading source code, making targeted edits to plugin_validator.py, and running tests. General-purpose agent has full tool access for code reading, editing, and test execution.

---

### #2: Validate and verify orchestrator-discipline plugin (P1)

**Rationale**: Fixes a broken P1 security/enforcement plugin. Plugin currently fails `claude plugin validate` — it is not functional as intended. Comprehensive plan pre-exists with all architectural decisions made, file paths specified, 5 hook test commands documented. T1/T2/T3 parallelizable; T4 convergence. One focused session can complete all 4 tasks.

**Dependencies**: None
**Prerequisites**: plan/tasks-4-validate-orchestrator-discipline.md (exists, 4 tasks ready)
**Complexity**: Medium — mix of plugin.json JSON editing, JavaScript hook modification, and SKILL.md frontmatter update

**Sub-agent**: `general-purpose`
**Why**: Multi-file task (plugin.json, .cjs hook, SKILL.md, rules/CLAUDE.md) requiring careful coordination. General-purpose agent can follow the plan file and execute all 4 tasks.

---

### #3: Add PR003/PR004 test coverage (P2)

**Rationale**: Smallest actionable item with zero blockers. Closes a concrete test coverage gap identified in code review. Existing PR001/PR002 test classes serve as copy-paste templates. Severity levels (info vs warning) already verified by groomer agent. No architectural decisions needed — pure test implementation.

**Dependencies**: None
**Prerequisites**: None (all patterns discovered by groomer agent)
**Complexity**: Small — Python test writing following established patterns

**Sub-agent**: `general-purpose`
**Why**: Pure Python test implementation with clear patterns. Groomer agent provided sample test class structures ready to adapt.

---

## Assignment Plan

### Item #1: plugin-validator UX and coverage gaps
**Agent**: general-purpose
**Instructions to sub-agent**:
- Work in `plugins/plugin-creator/scripts/plugin_validator.py`
- Follow order: Sub-issue 4 (remove lines 904-911 dead code + error strings at 759-760, 771-772) → Sub-issue 1 (fix reporter at lines 2928-2931 to count unique files, not validators) → Sub-issue 3 (note: FileType.HOOK_CONFIG already exists at lines 318-329; add HOOK_SCRIPT, add HookValidator class, add HK001+ codes) → Sub-issue 2 (add file_type check to DescriptionValidator, suppress SK005 for COMMAND files)
- Run `uv run pytest plugins/plugin-creator/tests/ -v` after each sub-issue
- Reference QA report at `plugins/plugin-creator/planning/plugin-validator-qa-report.md` for context on 15 pre-existing failures
- Reference plan at `plan/tasks-2-validator-ux-coverage.md`

### Item #2: Validate orchestrator-discipline plugin
**Agent**: general-purpose
**Instructions to sub-agent**:
- Read `plan/tasks-4-validate-orchestrator-discipline.md` in full before starting
- Execute T1 (fix plugin.json `rules` field), T2 (add directory detection to hook), T3 (add `user-invocable: true` to SKILL.md) in parallel if possible, then T4 (convergence validation)
- Run `claude plugin validate plugins/orchestrator-discipline/` after T1
- Run the 5 hook test commands specified in plan after T2
- Run `uv run prek run --files plugins/orchestrator-discipline/` after T3
- Document any pre-existing issues found in T4 output per plan line 800-801

### Item #3: Add PR003/PR004 test coverage
**Agent**: general-purpose
**Instructions to sub-agent**:
- Target file: `plugins/plugin-creator/tests/test_plugin_registration_validator.py`
- Read existing PR001 (lines 267-328) and PR002 (lines 470-543) test classes for patterns
- PR003 = info severity → check `result.info` list; emitted when metadata fields (repository, homepage, author) absent
- PR004 = warning severity → check `result.warnings` list; emitted when plugin.json repository differs from git remote URL
- Create `TestMissingMetadata` class for PR003 tests
- Create `TestRepositoryMismatch` class for PR004 tests — may need `monkeypatch` for git metadata
- Run `uv run pytest plugins/plugin-creator/tests/test_plugin_registration_validator.py -v`

---

## Backlog Corrections Needed

1. Update path for all ARL-related items: `skills/assessor/references/ARL/` → `skills/arl/references/`
2. Promote "Configurable Token Thresholds" from Idea to P2 (zero blockers, small scope)
3. Add note to "SAM: Replace validate-task-file.sh" — verify script is actually used before implementing replacement
4. Mark "Validate is-fast", "Validate agent-browser", "Validate carbonyl" as BLOCKED pending environment access
