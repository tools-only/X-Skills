# Refactoring Tasks: holistic-linting

**Plugin**: ./plugins/holistic-linting
**Design Spec**: [refactor-design-holistic-linting.md](./refactor-design-holistic-linting.md)
**Baseline Score**: 68/100
**Target Score**: 85+/100

---

## Task 1: Split holistic-linting SKILL.md into 3 focused skills

**Status**: ❌ NOT STARTED
**Dependencies**: None
**Priority**: 1
**Complexity**: High
**Agent**: plugin-creator:refactor-skill

**Target**: ./plugins/holistic-linting/skills/holistic-linting/SKILL.md
**Issue Type**: SKILL_SPLIT

**Description**:

Split the 1221-line monolithic SKILL.md into 3 focused skills:

1. **holistic-linting** (router/core) — Trimmed to ~430 lines containing: purpose (lines 1-48), linter detection and running formatters (lines 247-375), bundled resources (lines 993-1222). This remains the primary entry point. Remove the `<section ROLE_TYPE="orchestrator">` block (moved to orchestrator skill) and the linter-specific resolution workflows (moved to resolver skill). Add cross-references to the two new skills.

2. **holistic-linting-orchestrator** — New skill ~200 lines containing the orchestrator delegation workflow currently in lines 49-245. This is the content inside `<section ROLE_TYPE="orchestrator">` including: agent delegation patterns, step-by-step workflow, anti-patterns to avoid, workflow summary.

3. **holistic-linting-resolver** — New skill ~500 lines containing linter-specific resolution workflows currently in lines 376-991 (Ruff 380-513, Mypy 514-707, Pyright 709-957, Integration 958-991). Deduplicate the shared suppression gate sections (duplicated 3x), "Load python3-development skill" steps, and verification patterns into a common methodology section.

**Additional fixes during split**:
- Lines 1005-1011: Change `python holistic-linting/scripts/install-agents.py` to `uv run ./scripts/install-agents.py` for path consistency

**Acceptance Criteria**:
1. Main holistic-linting/SKILL.md is under 450 lines and contains only router/core content
2. holistic-linting-orchestrator/SKILL.md exists with ~200 lines of orchestrator-specific content
3. holistic-linting-resolver/SKILL.md exists with under 520 lines, deduplication applied
4. All three skills have valid YAML frontmatter with appropriate descriptions containing trigger keywords
5. Cross-references between the three skills use correct relative paths
6. No content is lost — all 1221 lines of original content are preserved across the three files (minus deduplication)
7. The install-agents.py path on former lines 1005-1011 uses `uv run ./scripts/install-agents.py`

**Required Inputs**:
- Design spec section: "Skill Split"
- Source files: ./plugins/holistic-linting/skills/holistic-linting/SKILL.md

**Expected Outputs**:
- ./plugins/holistic-linting/skills/holistic-linting/SKILL.md (trimmed)
- ./plugins/holistic-linting/skills/holistic-linting-orchestrator/SKILL.md (new)
- ./plugins/holistic-linting/skills/holistic-linting-resolver/SKILL.md (new)

**Can Parallelize With**: Task 2
**Reason**: Task 1 modifies skill files; Task 2 modifies an agent file — no shared files

**Verification Steps**:
1. Run `wc -l ./plugins/holistic-linting/skills/holistic-linting/SKILL.md` — under 450 lines
2. Run `wc -l ./plugins/holistic-linting/skills/holistic-linting-orchestrator/SKILL.md` — approximately 200 lines
3. Run `wc -l ./plugins/holistic-linting/skills/holistic-linting-resolver/SKILL.md` — under 520 lines
4. Verify frontmatter in all three files: `head -5 ./plugins/holistic-linting/skills/*/SKILL.md`
5. Grep for cross-references: `grep -r "holistic-linting-orchestrator\|holistic-linting-resolver" ./plugins/holistic-linting/skills/holistic-linting/SKILL.md`
6. Verify no suppression gate duplication: `grep -c "Suppression Gate" ./plugins/holistic-linting/skills/holistic-linting-resolver/SKILL.md` — should be 1 (shared section)
7. Verify install path fix: `grep "install-agents" ./plugins/holistic-linting/skills/holistic-linting/SKILL.md` — should show `uv run ./scripts/`

---

## Task 2: Optimize post-linting-architecture-reviewer agent

**Status**: ❌ NOT STARTED
**Dependencies**: None
**Priority**: 2
**Complexity**: Medium
**Agent**: subagent-refactorer

**Target**: ./plugins/holistic-linting/agents/post-linting-architecture-reviewer.md
**Issue Type**: AGENT_OPTIMIZE

**Description**:

Optimize the post-linting-architecture-reviewer agent:

1. **Change model from `haiku` to `inherit`** (line 4). Architectural review requires reliable reasoning about design patterns, type safety, and code organization. Haiku has documented ~50% hallucination rate for reasoning tasks (per project CLAUDE.md).

2. **Enhance description with additional trigger keywords**: Add "design principles", "type safety validation", "SOLID compliance", "code organization" to improve delegation matching accuracy.

3. **Fix markdown structure**: The code fence block starting at line 116 (````markdown) has structural issues. Lines 164-165 close the inner fences but line 182 has a stray closing ```. Clean up the fence nesting so the file renders correctly.

**Acceptance Criteria**:
1. Frontmatter `model` field is `inherit` (not `haiku`)
2. Frontmatter `description` contains keywords: "design principles", "type safety", "code organization"
3. No stray/unclosed markdown code fences in the file
4. File renders correctly as markdown (no broken fence blocks)
5. All existing review checklist items (lines 54-111) preserved unchanged

**Required Inputs**:
- Design spec section: "Agent Optimizations"
- Source files: ./plugins/holistic-linting/agents/post-linting-architecture-reviewer.md

**Expected Outputs**:
- ./plugins/holistic-linting/agents/post-linting-architecture-reviewer.md (modified)

**Can Parallelize With**: Task 1
**Reason**: Modifies agent file only; Task 1 modifies skill files — no shared files

**Verification Steps**:
1. Run `head -6 ./plugins/holistic-linting/agents/post-linting-architecture-reviewer.md` — verify model: inherit
2. Grep description for new keywords: `grep -i "design principles\|type safety\|code organization" ./plugins/holistic-linting/agents/post-linting-architecture-reviewer.md`
3. Validate markdown fence balance: Count opening and closing fences match
4. Verify review checklist preserved: `grep -c "\- \[ \]" ./plugins/holistic-linting/agents/post-linting-architecture-reviewer.md` — should be >=20

---

## Task 3: Update plugin.json with new skills and metadata

**Status**: ❌ NOT STARTED
**Dependencies**: Task 1
**Priority**: 2
**Complexity**: Low
**Agent**: claude-context-optimizer

**Target**: ./plugins/holistic-linting/.claude-plugin/plugin.json
**Issue Type**: STRUCTURE_FIX

**Description**:

Update plugin.json to:
1. Add the two new skill directories to the `skills` array
2. Add `keywords` field for marketplace discoverability
3. Bump version from 1.0.7 to 1.1.0 (minor version for new skills)

**Acceptance Criteria**:
1. `skills` array contains all 3 skill paths: `["./skills/holistic-linting", "./skills/holistic-linting-orchestrator", "./skills/holistic-linting-resolver"]`
2. `keywords` field present with relevant terms: linting, code-quality, ruff, mypy, bandit, python, formatting, type-checking
3. `version` bumped to `"1.1.0"`
4. File is valid JSON (parseable by `python3 -m json.tool`)

**Required Inputs**:
- Design spec section: "Documentation Improvements" items 2 and 3
- Source files: ./plugins/holistic-linting/.claude-plugin/plugin.json

**Expected Outputs**:
- ./plugins/holistic-linting/.claude-plugin/plugin.json (modified)

**Can Parallelize With**: None (depends on Task 1)
**Reason**: Needs Task 1 complete to reference correct skill directories

**Verification Steps**:
1. Run `python3 -m json.tool ./plugins/holistic-linting/.claude-plugin/plugin.json` — valid JSON
2. Verify skills array: `python3 -c "import json; d=json.load(open('./plugins/holistic-linting/.claude-plugin/plugin.json')); print(len(d['skills']))"` — should be 3
3. Verify keywords present: `grep keywords ./plugins/holistic-linting/.claude-plugin/plugin.json`
4. Verify version: `grep version ./plugins/holistic-linting/.claude-plugin/plugin.json` — should show 1.1.0

---

## Task V1: Validate plugin structure

**Status**: ❌ NOT STARTED
**Dependencies**: Task 1, Task 2, Task 3
**Priority**: 1
**Complexity**: Low
**Agent**: plugin-assessor

**Target**: ./plugins/holistic-linting/
**Issue Type**: STRUCTURE_FIX

**Acceptance Criteria**:
1. Plugin passes structural validation with no critical issues
2. All internal links resolve correctly (skills reference each other, agent references resolve)
3. No orphaned files (status maintained from baseline — 0 orphans)
4. Frontmatter validates against Claude Code schema for all 3 skills and 2 agents
5. Score improved from 68/100 baseline

**Required Inputs**:
- Source files: ./plugins/holistic-linting/ (entire plugin directory)

**Expected Outputs**:
- Validation report with updated score

**Can Parallelize With**: None
**Reason**: Depends on all modification tasks completing

**Verification Steps**:
1. Run plugin-assessor on ./plugins/holistic-linting/ — no critical issues
2. Verify all 3 skill SKILL.md files exist and have valid frontmatter
3. Verify both agent files have valid frontmatter
4. Check all markdown links resolve: grep for `](./` in all SKILL.md files and verify targets exist
5. Score shows improvement from 68/100

---

## Task V2: Update plugin documentation

**Status**: ❌ NOT STARTED
**Dependencies**: Task V1
**Priority**: 2
**Complexity**: Low
**Agent**: plugin-docs-writer

**Target**: ./plugins/holistic-linting/README.md
**Issue Type**: DOC_IMPROVE

**Acceptance Criteria**:
1. README.md reflects current plugin structure with 3 skills
2. All skills listed with descriptions matching their frontmatter
3. Both agents documented with their roles
4. The /lint command documented with usage examples
5. Installation instructions accurate

**Required Inputs**:
- Source files: ./plugins/holistic-linting/ (entire plugin directory)
- Validation results from Task V1

**Expected Outputs**:
- ./plugins/holistic-linting/README.md (updated)

**Can Parallelize With**: None
**Reason**: Depends on Task V1 (needs validated structure to document accurately)

**Verification Steps**:
1. README.md exists and has content
2. All 3 skill names appear in README: `grep -c "holistic-linting-orchestrator\|holistic-linting-resolver\|holistic-linting" ./plugins/holistic-linting/README.md`
3. Both agent names appear: `grep -c "linting-root-cause-resolver\|post-linting-architecture-reviewer" ./plugins/holistic-linting/README.md`
4. /lint command documented: `grep "/lint" ./plugins/holistic-linting/README.md`

---

## Dependency Summary

```text
Phase 1 (Parallel):
  Task 1: Split SKILL.md ──────────┐
  Task 2: Optimize agent ──────────┤
                                   │
Phase 2 (Sequential):              │
  Task 3: Update plugin.json ◄─────┘ (depends on Task 1)
                                   │
Phase 3 (Sequential):              │
  Task V1: Validate structure ◄────┘ (depends on Tasks 1, 2, 3)
                                   │
Phase 4 (Sequential):              │
  Task V2: Update README ◄────────┘ (depends on Task V1)
```

## Parallelization Groups

- **Group A** (can launch simultaneously): Task 1, Task 2 — No shared files
- **Group B** (sequential after Group A): Task 3 — Depends on Task 1
- **Group C** (sequential after all modifications): Task V1 — Depends on Tasks 1, 2, 3
- **Group D** (sequential after validation): Task V2 — Depends on Task V1

## Success Metrics

| Metric | Baseline | Target |
|--------|----------|--------|
| Overall plugin score | 68/100 | 85+/100 |
| Largest skill file | 1221 lines | <520 lines |
| Number of skills | 1 | 3 |
| Agent model appropriateness | haiku for reasoning | inherit |
| Critical issues | 0 (corrected from assessment) | 0 |
| Orphaned files | 0 | 0 |

---

## Context Manifest

### Domain Boundaries in SKILL.md

The source SKILL.md (1221 lines) contains 4 distinct domains that will be split across 3 skills:

| Domain | Line Range | Target Skill | Content Summary |
|--------|-----------|--------------|-----------------|
| Purpose, Overview, Linter Detection | 1-48 | holistic-linting (router) | Purpose statement, when skill applies, orchestrator vs sub-agent roles, linter detection table |
| Orchestrator Delegation Workflow | 49-245 | holistic-linting-orchestrator | Inside `<section ROLE_TYPE="orchestrator">` tag. Complete delegation workflow, agent invocation patterns, anti-patterns, report reading, multi-file handling |
| Running Formatters and Linters | 247-375 | holistic-linting (router) | How to run formatters/linters, git hook tool detection, scoped operations, language-specific commands |
| Linter-Specific Resolution Workflows | 376-991 | holistic-linting-resolver | Ruff workflow (380-513), Mypy workflow (514-707), Pyright workflow (709-957), Integration notes (958-991). Contains suppression gate (duplicated 3x), load python3-development step (duplicated 3x), architectural context check (similar 3x) |
| Bundled Resources, Examples, Best Practices | 993-1222 | holistic-linting (router) | Agent installation instructions, rules knowledge base documentation, /lint command reference, hooks integration, examples, troubleshooting |

**Deduplication Opportunity**: Lines 376-991 contain ~36 lines of duplicated suppression gate sections (12 lines × 3 linters) and shared methodology steps. Extract common resolution pattern into shared section to reduce from ~616 lines to target ~500 lines.

**Path Fix**: Lines 1005-1011 use `python holistic-linting/scripts/install-agents.py` which is relative to skills directory parent. Will be changed to `uv run ./scripts/install-agents.py` for consistency with other script references.

### Cross-Reference Map

| Source File | References To | Reference Type |
|-------------|--------------|----------------|
| SKILL.md line 997 | `./agents/linting-root-cause-resolver.md` | Markdown link to bundled agent |
| SKILL.md line 1027 | `./references/rules/ruff/index.md` | Markdown link to rules KB index |
| SKILL.md line 1050 | `./references/rules/mypy/index.md` | Markdown link to rules KB index |
| SKILL.md line 1074 | `./references/rules/bandit/index.md` | Markdown link to rules KB index |
| SKILL.md line 1100 | `./scripts/` | Directory reference for bundled scripts |
| SKILL.md line 1118 | `/.claude/commands/lint.md` | Absolute path to command (works from plugin context) |
| SKILL.md lines 453, 593, 821 | `Skill(command: "python3-development")` | Skill activation in all 3 linter workflows |
| linting-root-cause-resolver.md line 17 | `Skill(command: "holistic-linting")` | Agent loads main skill |
| linting-root-cause-resolver.md line 22 | `Skill(command: "python3-development")` | Agent loads python skill |
| post-linting-architecture-reviewer.md line 29 | `.claude/reports/linting-investigation-*.md` | Resolution artifact path |
| post-linting-architecture-reviewer.md line 30 | `.claude/reports/linting-resolution-*.md` | Resolution artifact path |
| post-linting-architecture-reviewer.md line 31 | `.claude/artifacts/linting-artifacts-*.json` | Resolution artifact path |
| lint.md line 62 | `Grep(pattern="^## LINTERS", path="CLAUDE.md")` | Command reads project linting config |
| lint.md line 114 | `Task(subagent_type="linting-root-cause-resolver")` | Command delegates to agent |

**After Split Impact**: The main `holistic-linting` skill remains the entry point. Agent references `Skill(command: "holistic-linting")` will continue to work as the router skill will reference or incorporate resolver content. No agent file changes needed for skill references.

### External Dependencies

External files that reference holistic-linting (refactoring will NOT break these):

| External File | References | Impact of Refactoring |
|--------------|-----------|----------------------|
| `.claude-plugin/marketplace.json` line 45-46 | Plugin entry: `"name": "holistic-linting", "source": "./plugins/holistic-linting"` | Plugin name unchanged, no impact |
| `plugins/python3-development/skills/python3-development/SKILL.md` line 828 | "per holistic-linting skill" | Generic reference to skill name, no impact |
| `plugins/python3-development/skills/python3-development/SKILL.md` line 984 | Path to detection script (informational) | Path remains valid, no impact |
| `plugins/agent-orchestration/skills/agent-orchestration/SKILL.md` lines 223, 258, 692, 927 | "per the holistic-linting skill" | Generic reference to skill name, no impact |
| `plugins/agent-orchestration/skills/agent-orchestration/SKILL.md` lines 703, 938 | "Activate holistic-linting skill" | Main skill name unchanged, no impact |
| `plugins/python3-development/agents/code-reviewer.md` | References linting workflows | Generic reference, no impact |
| `plugins/verification-gate/skills/verification-gate/SKILL.md` | References holistic linting concept | Generic reference, no impact |
| `README.md` | Lists plugin in repository overview | Plugin name unchanged, no impact |

**Assessment**: No external file changes required. All references are either to the plugin name (unchanged) or generic mentions of "holistic-linting skill" which remain valid as the main skill persists as the entry point.

### Reference Index Structure

The rules knowledge base uses a three-tier progressive disclosure hierarchy:

**Tier 1: SKILL.md** (lines 1021-1097) — Documents existence of rules knowledge base with summary tables:

- Ruff: 933 rules, 19 families documented
- Mypy: Error codes organized by category
- Bandit: 65+ security checks organized by vulnerability type

**Tier 2: Index Files** (3 files) — Provide rule family/category organization and navigation:

- `./references/rules/ruff/index.md` — Lists 19 rule families (pycodestyle, pyflakes, bugbear, etc.) with links to family-specific files
- `./references/rules/mypy/index.md` — Organizes error codes by category (attribute access, function calls, etc.) with links to error documentation
- `./references/rules/bandit/index.md` — Organizes security checks by vulnerability category (credentials, crypto, injection, etc.)

**Tier 3: Detailed Rule Files** (52 files total in `./references/rules/{ruff,mypy,bandit}/`) — Individual rule documentation with examples, principles, and resolution patterns

**Link Patterns**: All links use relative markdown paths starting with `./` from the referencing file's location. SKILL.md → `./references/rules/ruff/index.md`, index files → `./pycodestyle-errors-warnings.md`

**After Split**: The router skill (holistic-linting) will retain documentation of the rules knowledge base (lines 1021-1097). The resolver skill (holistic-linting-resolver) will reference the rules KB but does NOT need to duplicate index links — agents activate the main skill which provides navigation to the KB.

### Additional Reference Files

Beyond the rules knowledge base, the plugin contains:

**Mypy Cached Documentation** (6 `.rst` files in `./references/mypy-docs/`):

- `error_code_list.rst` and `error_code_list2.rst` — Referenced in SKILL.md lines 526-529 as locally-cached mypy error code documentation
- Used by Mypy Resolution Workflow for offline rule lookup
- Paths: `./references/mypy-docs/error_code_list.rst` (relative from SKILL.md)

**Online Resource Access Guide**:

- `./references/accessing_online_resources.md` — Referenced in SKILL.md line 749 for MCP tool usage guidance when researching Pyright rules
- Provides instructions for using `mcp__Ref__ref_search_documentation`, `mcp__exa__get_code_context_exa`, WebFetch tools

**Scripts**:

- `./scripts/detect-hook-tool.py` — Git hook tool detection (prek vs pre-commit)
- `./scripts/install-agents.py` — Agent installation to user/project scope
- Referenced in SKILL.md lines 302-310 (detection) and 1005-1011 (installation)

All reference files use relative paths from SKILL.md location and remain unchanged by the split.

### Key Observations

**Architectural Pattern**: The current SKILL.md functions as a "kitchen sink" containing:

1. Router logic (what skill does, when to use, linter detection)
2. Role-specific workflows (orchestrator delegation inside XML section tag)
3. Execution procedures (linter-specific resolution workflows for sub-agents)
4. Resource documentation (bundled agents, rules KB, scripts)

**Split Rationale**: The XML section tag on line 49 (`<section ROLE_TYPE="orchestrator">`) explicitly delineates orchestrator-only content, suggesting original intent to separate concerns. The refactoring actualizes this separation into discrete skills.

**Deduplication Targets**:

- Suppression Gate section: Lines 458-470 (Ruff), 598-610 (Mypy), 825-837 (Pyright) — Identical text duplicated 3 times
- "Load python3-development skill" step: Lines 449-456 (Ruff), 589-596 (Mypy), 816-823 (Pyright) — Near-identical with slight wording variations
- Verification step: Lines 483-488 (Ruff), 668-674 (Mypy), 911-919 (Pyright) — Similar pattern with tool name variation

Extract common methodology once, then linter-specific sections reference it with tool-specific variations.

**Agent Model Issue**: `post-linting-architecture-reviewer.md` line 4 sets `model: haiku`. Architectural review requires reasoning about design patterns, code organization, type safety, SOLID principles — tasks where Haiku has documented ~50% hallucination rate per project CLAUDE.md. Model should be `inherit` to use orchestrator's reasoning-capable model.

**Markdown Structure Issue**: `post-linting-architecture-reviewer.md` lines 116-182 have fence nesting problems. Opening `````markdown` on line 116, but closing structure is malformed with stray backticks on lines 165, 182. This prevents proper markdown rendering.

**No Orphaned Files**: All 52 reference files in `./references/rules/` are indexed through the three-tier hierarchy. The 6 mypy-docs RST files are explicitly referenced in SKILL.md. The accessing_online_resources.md is linked inline. Scripts are documented in bundled resources section. Assessment confirmed: 0 orphaned files.

**Version Bump Required**: plugin.json version is `1.0.7`. Adding two new skills constitutes a minor feature addition. Version should bump to `1.1.0` per repository version bumping rules.

**Cross-Plugin Integration Points**:

- python3-development skill: Activated in all 3 linter resolution workflows (lines 453, 593, 821) to ensure Python 3.11+ standards
- agent-orchestration skill: References holistic-linting workflows in delegation templates
- verification-gate skill: May reference linting validation in completion gates

These integration points use skill activation syntax (`Skill(command: "holistic-linting")`) which continues to work as the main skill name is unchanged.
