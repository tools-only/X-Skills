# Feature Context: Enhance skill-research-process for Complete CLI Tool Skill Production

## Document Metadata

- **Generated**: 2026-02-27
- **Input Type**: assessment_gaps (skill comparison against reference benchmark)
- **Source**: Session observation comparing skill-research-process output against uv skill; Issue #197
- **Status**: DISCOVERY_COMPLETE
- **Priority**: P1

---

## Original Request

The skill-research-process skill has sound research orchestration but lacks output specification for producing complete CLI tool skills. Five gaps identified via assessment against the uv skill:

1. No local directory input -- passing a path like `.claude/worktrees/ty/docs/` is treated as a tool name, triggering web searches instead of reading local docs
2. No CLI reference file templates -- no structural anchor ensuring standard reference types (`cli_reference.md`, `configuration.md`, `migration-guide.md`, `quick-reference.md`, `troubleshooting.md`)
3. No `assets/` directory production
4. No sync/release-tracking script production
5. Output uses `references/{category}/index.md` (subdirectory) instead of flat `references/{slug}.md` layout

---

## Core Intent Analysis

### WHO (Target Users)

- Skill authors using `/skill-research-process` to build new skills for CLI tools, frameworks, and libraries
- Orchestrator agents delegating research work to parallel sub-agents
- End users consuming the produced skills (who benefit from consistent structure)

### WHAT (Desired Outcome)

- The skill-research-process accepts local documentation directories as input alongside tool names
- Research output conforms to the proven reference layout observed in production-quality skills (uv, clang-format, gitlab-skill)
- Produced skills include `assets/` directory with templates, `scripts/` with sync scripts, and flat `references/*.md` files
- A single invocation of `/skill-research-process` produces a skill structure that passes `plugin_validator.py` without post-hoc restructuring

### WHEN (Trigger Conditions)

- User invokes `/skill-research-process <path-to-local-docs>` (local directory input)
- User invokes `/skill-research-process <tool-name>` (web-based research, existing flow)
- Research agents complete Stage 2 and produce reference files (output layout enforcement)
- Integration agent runs Stage 3 (assets/ and scripts/ production)

### WHY (Problem Being Solved)

- Local docs (cloned repos, worktrees, downloaded archives) are a high-fidelity source that the current skill ignores, forcing unnecessary web searches
- Without structural templates, each research run produces ad-hoc reference layouts that require manual restructuring to match repository conventions
- Missing `assets/` and `scripts/` means the produced skill is incomplete compared to the repository's quality benchmark (uv skill)
- Subdirectory layout (`references/{category}/index.md`) diverges from the flat layout (`references/{slug}.md`) used by every production skill in the repository

---

## Codebase Research

### Similar Patterns Found

#### Pattern 1: Local path input handling — external-pattern-integrator

- **Location**: `.claude/skills/external-pattern-integrator/SKILL.md:68-72`
- **Code**:

  ```text
  **If URL**: Use WebFetch or curl to download to `/tmp/external-pattern-{slug}.md`
  **If local file**: Read directly
  ```

- **Relevance**: Demonstrates argument-level branching between URL and local file input within a skill
- **Reusable**: The `argument-hint: <url-or-file> [url-or-file...]` pattern and the if-URL/if-local branching logic

#### Pattern 2: Flat references/ layout — uv skill (benchmark)

- **Location**: `plugins/python3-development/skills/uv/references/`
- **Files**:
  - `cli_reference.md`
  - `configuration.md`
  - `migration-guide.md`
  - `quick-reference.md`
  - `troubleshooting.md`
- **Relevance**: Production-quality CLI tool skill uses flat `references/{slug}.md` -- no subdirectories, no `index.md` files
- **Reusable**: File naming convention and standard reference types for CLI tools

#### Pattern 3: Flat references/ layout — brainstorming-skill

- **Location**: `plugins/brainstorming-skill/skills/brainstorming-skill/SKILL.md:73-108`
- **Code**: References 14+ files all as `./references/{slug}.md` (flat layout)
- **Relevance**: Non-CLI skill also uses flat layout, confirming this is the repository-wide convention
- **Reusable**: Confirms the pattern is not CLI-specific

#### Pattern 4: Flat references/ layout — agent-browser skill

- **Location**: `.claude/skills/agent-browser/SKILL.md:393-399`
- **Code**: `references/commands.md`, `references/snapshot-refs.md`, `references/session-management.md`, etc.
- **Relevance**: Another production skill using flat references
- **Reusable**: Confirms flat layout as universal standard

#### Pattern 5: assets/ directory with categorized subdirectories — uv skill

- **Location**: `plugins/python3-development/skills/uv/assets/`
- **Files**:
  - `pyproject_templates/basic.toml`, `advanced.toml`, `gitlab.toml`
  - `script_examples/data_analysis.py`
  - `docker_examples/Dockerfile.simple`, `Dockerfile.multi-stage`
  - `github_actions/ci.yml`
- **Relevance**: Assets organized by purpose subdirectory, providing copy-paste templates for users
- **Reusable**: Asset category structure (templates, examples, configs)

#### Pattern 6: assets/ directory with categorized subdirectories — clang-format skill

- **Location**: `plugins/clang-format/skills/clang-format/assets/`
- **Files**:
  - `configs/` (7 `.clang-format` configuration files)
  - `integrations/` (emacs, vim, pre-commit integrations)
- **Relevance**: Another production skill with meaningful assets
- **Reusable**: Config templates and integration examples as asset categories

#### Pattern 7: Sync/release-tracking script — uv skill

- **Location**: `plugins/python3-development/skills/uv/scripts/sync_uv_releases.py`
- **Functionality**: Fetches GitHub releases, categorizes changes (breaking, features, deprecations), updates SKILL.md Version Information section, includes cooldown with lock file
- **Relevance**: Directly demonstrates the sync script pattern missing from skill-research-process output
- **Reusable**: Script structure (PEP 723 metadata, typer CLI, cooldown logic, lock file, section replacement in SKILL.md)

#### Pattern 8: Sync script — gitlab-skill

- **Location**: `plugins/gitlab-skill/skills/gitlab-skill/scripts/sync_gitlab_docs.py`
- **Functionality**: Downloads documentation archive, extracts, grooms markdown (link transforms, Hugo shortcode removal), generates file tree index, atomic replacement
- **Relevance**: Alternative sync pattern (archive download vs API) for documentation-heavy skills
- **Reusable**: Archive-based sync for tools with static documentation sites

#### Pattern 9: Add-doc-updater orchestration — skill-creator plugin

- **Location**: `plugins/plugin-creator/skills/add-doc-updater/SKILL.md`
- **Functionality**: 5-phase workflow to add documentation sync pipeline to any skill; collects 6 template variables; delegates to `@python-cli-architect` agent
- **Relevance**: Provides existing infrastructure for adding sync scripts to skills -- skill-research-process could invoke this as a Stage 3 step
- **Reusable**: The `/add-doc-updater` skill itself; no need to reinvent sync script generation

#### Pattern 10: Skill scaffolding with all three directories — init_skill.py

- **Location**: `plugins/plugin-creator/skills/skill-creator/scripts/init_skill.py:259-307`
- **Code**: Creates `scripts/example.py`, `references/api_reference.md`, `assets/example_asset.txt`
- **Relevance**: The scaffolding script already creates all three directories; skill-research-process should populate them instead of leaving assets/ and scripts/ empty
- **Reusable**: Directory creation pattern and placeholder structure

### Existing Infrastructure

| Component | Path | Role |
|-----------|------|------|
| skill-research-process | `.claude/skills/skill-research-process/SKILL.md` | Current skill (target of enhancement) |
| Agent prompts | `.claude/skills/skill-research-process/references/agent-prompts.md` | Research/categorization/integration agent templates |
| Gaps analysis | `.claude/skills/skill-research-process/references/gaps-analysis.md` | Known gaps (9 items, but none cover the 5 gaps in this request) |
| MCP tools guide | `.claude/skills/skill-research-process/references/mcp-tools.md` | Tool selection reference |
| init_skill.py | `plugins/plugin-creator/skills/skill-creator/scripts/init_skill.py` | Skill scaffolding |
| add-doc-updater | `plugins/plugin-creator/skills/add-doc-updater/SKILL.md` | Sync script generator |
| plugin_validator.py | `plugins/plugin-creator/scripts/plugin_validator.py` | Structural validation |

### Code References

- `.claude/skills/skill-research-process/SKILL.md:2` -- `argument-hint: <tool-or-library-name>` (current: tool name only, no path)
- `.claude/skills/skill-research-process/SKILL.md:21` -- `Stage 2: Research → Parallel agents populate references/{category}/` (subdirectory layout)
- `.claude/skills/skill-research-process/SKILL.md:75` -- `Each agent outputs to ./references/{category}/` (subdirectory layout)
- `.claude/skills/skill-research-process/references/agent-prompts.md:64` -- `Create files in: ./{skill-name}/references/{category}/` (subdirectory layout in agent prompt)
- `.claude/skills/skill-research-process/references/agent-prompts.md:65` -- `Create index.md (lowercase) in that directory` (index.md pattern)
- `.claude/skills/skill-research-process/SKILL.md:112` -- `Update ./SKILL.md with links to each category's index.md` (integration references index.md)
- `plugins/python3-development/skills/uv/references/` -- 5 flat `.md` files, zero subdirectories (benchmark)
- `plugins/python3-development/skills/uv/assets/` -- 4 subdirectories with 7 template files (benchmark)
- `plugins/python3-development/skills/uv/scripts/sync_uv_releases.py` -- 587-line sync script (benchmark)

---

## Use Scenarios

### Scenario 1: Building a skill from cloned repository docs

**Actor**: Skill author who has cloned `astral-sh/ty` and wants to create a `/ty` skill
**Trigger**: `/skill-research-process .claude/worktrees/ty/docs/`
**Goal**: Research agents read local markdown files instead of web-searching "ty"
**Expected Outcome**: Categorization agent scans the local docs directory structure, creates categories from directory layout, research agents read local files with `Read` tool
**Current State**: The argument `.claude/worktrees/ty/docs/` is interpreted as a tool name. Categorization agent web-searches for "ty docs" -- returning unrelated results or nothing useful.

### Scenario 2: Building a CLI tool skill with complete output

**Actor**: Skill author building a `/kubectl` skill
**Trigger**: `/skill-research-process kubectl`
**Goal**: Produce a skill with flat references (`cli_reference.md`, `configuration.md`, `troubleshooting.md`), assets (YAML templates), and a sync script
**Expected Outcome**: Stage 2 produces `references/cli_reference.md`, `references/configuration.md`, `references/quick-reference.md`, `references/troubleshooting.md`; Stage 3 produces `assets/` with example configs and `scripts/sync_kubectl_releases.py`
**Current State**: Stage 2 produces `references/installation/index.md`, `references/commands/index.md`, etc. (subdirectory layout). No assets/ or scripts/ are produced. Manual restructuring required post-research.

### Scenario 3: Hybrid input -- local docs plus web enrichment

**Actor**: Skill author with partial local documentation wanting to supplement from web
**Trigger**: `/skill-research-process --local .claude/worktrees/ty/docs/ --name ty`
**Goal**: Research agents prioritize local files, fall back to web for topics not covered locally
**Expected Outcome**: Categorization scans local docs first, identifies gaps, research agents use local files where available and web sources where not
**Current State**: No mechanism to combine local and web sources. Must choose one or the other (and local is not supported at all).

---

## Gap Analysis

### Identified Gaps

| # | Category | Gap Description | Impact | Benchmark Reference |
|---|----------|-----------------|--------|---------------------|
| 1 | Input | No local directory path input -- argument treated as tool name | High-fidelity local docs cannot be used; forces unnecessary web searches with lower accuracy | `external-pattern-integrator` SKILL.md:68-72 handles URL-vs-file branching |
| 2 | Output Structure | No CLI reference file templates -- no guarantee of standard files like `cli_reference.md`, `configuration.md` | Each research run produces different file names; manual post-hoc rename required | uv skill: 5 named reference files at `references/*.md` |
| 3 | Output Structure | No `assets/` directory production | Produced skill lacks templates, example configs, CI workflows that users copy-paste | uv skill: 7 asset files across 4 subdirectories |
| 4 | Output Structure | No sync/release-tracking script production | Skill documentation becomes stale; no automated refresh mechanism | uv skill: `scripts/sync_uv_releases.py` (587 lines); existing `/add-doc-updater` skill |
| 5 | Output Structure | Subdirectory layout (`references/{category}/index.md`) instead of flat (`references/{slug}.md`) | Diverges from every production skill in the repository; requires restructuring before merge | All production skills: uv, clang-format, agent-browser, brainstorming-skill use flat layout |
| 6 | Agent Prompts | Research agent prompt (agent-prompts.md:64-65) hardcodes subdirectory + index.md pattern | Agents produce wrong layout; orchestrator cannot simply change SKILL.md without also updating prompts | agent-prompts.md lines 64-68 |
| 7 | Agent Prompts | Categorization agent prompt has no local-docs scanning path | Even if input parsing is fixed, the categorization agent does not know how to scan a directory | agent-prompts.md lines 10-42 |
| 8 | Integration | Stage 3 integration agent prompt (agent-prompts.md:134) references `references/{category}/index.md` | Integration step wires wrong paths into SKILL.md | agent-prompts.md line 134 |

---

## Questions Requiring Resolution

### Q1: How should local path vs. tool name be distinguished in the argument?

- **Category**: Input
- **Gap**: #1
- **Question**: Should the skill detect paths by checking if the argument is an existing directory, or require an explicit flag like `--local`?
- **Options**:
  - A) Auto-detect: If `$ARGUMENTS` resolves to an existing directory, treat as local docs path; otherwise treat as tool name
  - B) Explicit flag: `--local <path>` for local docs, bare argument for tool name
  - C) Both: Auto-detect with optional `--local` override for disambiguation
- **Why It Matters**: Auto-detect is simpler but may misfire on tool names that coincidentally match directory names. Explicit flag is unambiguous but adds syntax.
- **Recommendation**: Option A (auto-detect) with fallback documentation. The `external-pattern-integrator` uses auto-detect (URL vs file) and this has worked reliably. A simple `Path($ARGUMENTS).is_dir()` check in the skill body or a pre-step Bash command suffices.
- **Resolution**: _pending_

### Q2: Should the standard reference file set be hardcoded or configurable per tool type?

- **Category**: Output Structure
- **Gap**: #2
- **Question**: Should the skill always produce the same 5-6 reference files, or should the categorization agent determine which files to produce?
- **Options**:
  - A) Hardcoded template: Always produce `cli_reference.md`, `configuration.md`, `migration-guide.md`, `quick-reference.md`, `troubleshooting.md`
  - B) Configurable: Categorization agent proposes files; quality gate verifies minimum set
  - C) Hybrid: Required minimum set (cli_reference, configuration, troubleshooting) + agent-proposed additions
- **Why It Matters**: Hardcoded ensures consistency but may not fit all tools (e.g., a library may not have CLI reference). Configurable preserves flexibility but risks inconsistency.
- **Recommendation**: Option C. Require a minimum set for CLI tools, allow the categorization agent to add domain-specific files. The quality gate at Stage 1 should verify the minimum set is present.
- **Resolution**: _pending_

### Q3: Should assets/ production be part of skill-research-process or delegated to a follow-up skill?

- **Category**: Output Structure
- **Gap**: #3
- **Question**: Should research agents produce assets, or should assets be created in a separate post-research step?
- **Options**:
  - A) Inline: Research agents produce assets alongside references during Stage 2
  - B) Post-research: A dedicated Stage 3.5 produces assets from the researched content
  - C) Delegated: Invoke `/add-doc-updater` or similar after research completes
- **Why It Matters**: Research agents focus on documentation extraction; asset creation requires different judgment (what templates are useful, what configs to include). Mixing concerns may reduce quality.
- **Recommendation**: Option B. Add a post-research asset generation step in Stage 3 (Integration) where the integration agent identifies template-worthy content from references and creates assets. This keeps research agents focused.
- **Resolution**: _pending_

### Q4: Should sync script production be integrated or delegated to /add-doc-updater?

- **Category**: Output Structure
- **Gap**: #4
- **Question**: Should skill-research-process generate sync scripts directly, or invoke the existing `/add-doc-updater` skill?
- **Options**:
  - A) Integrated: Skill-research-process generates the sync script as part of Stage 3
  - B) Delegated: Stage 3 invokes `/add-doc-updater <produced-skill-path>` as a final step
  - C) Optional: Stage 3 asks user if sync script is needed (not all tools have release APIs)
- **Why It Matters**: `/add-doc-updater` already has a 5-phase workflow with quality gates for sync script creation. Reimplementing this inside skill-research-process duplicates effort and misses existing quality gates.
- **Recommendation**: Option B. The `/add-doc-updater` skill at `plugins/plugin-creator/skills/add-doc-updater/SKILL.md` already handles this with a complete 5-phase pipeline. Stage 3 should note "invoke `/add-doc-updater` as follow-up" rather than reimplementing.
- **Resolution**: _pending_

### Q5: Should the flat layout migration update the existing gaps-analysis.md or replace it?

- **Category**: Documentation
- **Gap**: #5, #6, #7, #8
- **Question**: The existing `references/gaps-analysis.md` documents 9 gaps. This feature adds 5 more. Should they be merged?
- **Options**:
  - A) Merge: Add the 5 new gaps to the existing gaps-analysis.md
  - B) Replace: Supersede with a new version that includes all gaps
  - C) Separate: Keep gaps-analysis.md for original gaps; add a new file for output-structure gaps
- **Why It Matters**: A single source of truth for gaps is easier to track. However, the original 9 gaps are about process quality (verification, citations, hallucination checks) while the new 5 are about output structure. Different concerns.
- **Recommendation**: Option A. Merge into a single gaps-analysis.md with two sections: "Process Quality Gaps" (original 9) and "Output Structure Gaps" (new 5). This keeps one file as the canonical gap tracker.
- **Resolution**: _pending_

---

## Goals (Pending Resolution)

_These goals will be finalized after questions are resolved._

1. **Update argument handling** in `SKILL.md` frontmatter and Stage 1 to detect and branch on local directory paths vs. tool names
2. **Update categorization agent prompt** in `references/agent-prompts.md` to support local directory scanning with `Read`/`Glob` tools instead of web search
3. **Replace subdirectory layout** (`references/{category}/index.md`) with flat layout (`references/{slug}.md`) across all three agent prompts and SKILL.md Stage 2/Stage 3 references
4. **Add CLI reference file template set** as a quality gate in Stage 1 -- categorization agent must include minimum reference types for CLI tools
5. **Add assets/ production step** in Stage 3 (Integration) where integration agent extracts template-worthy content into `assets/` subdirectories
6. **Add sync script delegation** in Stage 3 -- note to invoke `/add-doc-updater` as follow-up step for tools with release APIs or updatable documentation
7. **Update gaps-analysis.md** to include the 5 output-structure gaps alongside the 9 existing process-quality gaps
8. **Validate** produced skill structure passes `plugin_validator.py` without warnings related to layout

---

## Next Steps

After questions are resolved:

1. Update "Resolution" fields in Questions section
2. Finalize Goals section with accepted options
3. Proceed to architecture design (modifications to SKILL.md, agent-prompts.md, and potentially gaps-analysis.md)
4. Create task decomposition with file-level edit targets
5. Execute implementation
