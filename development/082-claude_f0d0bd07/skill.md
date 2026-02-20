# Plugin-Creator Plugin - AI-Facing Documentation

Complete plugin development toolkit for creating, refactoring, and validating Claude Code plugins, agents, skills, and commands.

---

## Plugin Identity

**Name:** `plugin-creator`
**Version:** 2.3.0
**Purpose:** Provides end-to-end capabilities for plugin development lifecycle - creation, validation, refactoring, and quality assurance.

**Core Capabilities:**

- ✅ Create new plugins (via script)
- ✅ Create new agents (via skill)
- ❌ Create new skills (GAP - no dedicated tool)
- ❌ Create new commands (GAP - no dedicated tool)
- ✅ Refactor plugins and skills
- ✅ Validate frontmatter, plugins, and structure
- ✅ Quality analysis and improvement

---

## CRITICAL: Skill Name Field Bug (Plugin Skills Only)

**Discovery Date:** 2026-01-29
**Status:** Confirmed Claude Code v2.1.23 bug
**Scope:** **Only affects plugin skills** - does NOT affect `.claude/skills/` directories

### The Bug

**Plugin skills** with an explicit `name:` field in their frontmatter **DO NOT appear as slash commands**, even when `user-invocable: true` is set.

**Important:** This bug only affects skills distributed via plugins. Skills in `.claude/skills/` (personal/project) work fine with or without the `name:` field.

**Symptoms (Plugin Skills Only):**

- Skill listed in `skills` array of plugin.json
- Skill has `user-invocable: true` in frontmatter
- Skill has `name: skill-name` in frontmatter
- **Result:** Skill does NOT appear in `/plugin-name:skill-name` autocomplete menu

**Not Affected:**

- Skills in `~/.claude/skills/` (user-level)
- Skills in `.claude/skills/` (project-level)
- These work fine with `name:` field present

**Workaround:**

- Remove the `name:` field entirely from SKILL.md frontmatter
- Claude Code will use the directory name instead
- **Result:** Skill WILL appear as `/plugin-name:skill-name`

### Evidence

**Test performed:**

1. plugin-creator had 14 skills, only `/count-lines` (a command) appeared
2. Removed `name:` field from 3 skills (skill-creator, agent-creator, assessor)
3. Those 3 skills immediately appeared as slash commands
4. Removed `name:` from all remaining skills
5. All skills now appear in autocomplete

**Before fix:**

```yaml
---
name: skill-creator  # <-- This PREVENTS slash command registration
description: Guide for creating effective skills
user-invocable: true
---
```

Result: `/plugin-creator:skill-creator` NOT in menu

**After fix:**

```yaml
---
description: Guide for creating effective skills
user-invocable: true
---
```

Result: `/plugin-creator:skill-creator` appears in menu ✓

### Validator Auto-Fix

The `plugin_validator.py` script now automatically removes `name:` fields from SKILL.md files with this explanation:

```
Removed 'name' field (Claude Code bug: skills with 'name' field don't appear as slash commands)
```

### Official Documentation Contradiction

The official Claude Code documentation states:

> `name`: Display name for the skill. If omitted, uses the directory name.

This implies `name:` is optional and should work either way. However, **having the field prevents slash command registration entirely**. This is a bug in Claude Code, not intended behavior.

### Recommendation

**For plugin skills:** DO NOT use `name:` field in SKILL.md files until this bug is fixed in Claude Code.

**For `.claude/skills/` (personal/project skills):** The `name:` field works correctly and can be used if desired (though it's optional since directory name is used by default).

---

## When to Use This Plugin

The model MUST use this plugin when:

- User asks to "create a new plugin"
- User asks to "create a new agent"
- User asks to "validate frontmatter" in SKILL.md or agent files
- User asks to "refactor a plugin" or "split a skill"
- User asks to "check skill size" or "check skill complexity"
- User needs to validate plugin.json or plugin structure
- User wants to fix tool formatting issues in frontmatter

---

## Component Inventory

### Skills (13)

| Skill                            | User-Invocable | Purpose                                                                                                     | Verified |
| -------------------------------- | -------------- | ----------------------------------------------------------------------------------------------------------- | -------- |
| `/plugin-creator`                | Yes            | Orchestrates plugin creation by delegating to specialist agents                                             | ✅ Yes   |
| `/agent-creator`                 | Yes            | Create agents from scratch or templates, handles scope (project/user/plugin), updates plugin.json if needed | ✅ Yes   |
| `/claude-skills-overview-2026`   | Yes            | Complete reference for Claude Code skills system (January 2026)                                             | ✅ Yes   |
| `/claude-plugins-reference-2026` | Yes            | Complete reference for Claude Code plugins system (January 2026)                                            | ✅ Yes   |
| `/claude-hooks-reference-2026`   | Yes            | Complete reference for Claude Code hooks system (January 2026)                                              | ✅ Yes   |
| `/add-doc-updater`               | Yes            | Orchestrate adding doc sync pipeline to skills wrapping external documentation (APIs, frameworks, CLIs). Creates Python script with cooldown, validation | ✅ Yes   |
| `/assessor`                      | Yes            | Analyze plugin structure and create refactoring task files                                                  | ✅ Yes   |
| `/ensure-complete`               | Yes            | Validate refactoring completeness and create follow-up tasks                                                | ✅ Yes   |
| `/feature-discovery`             | No             | Research features and identify gaps (delegated)                                                             | ✅ Yes   |
| `/implement-refactor`            | Yes            | Execute refactoring tasks with parallel orchestration                                                       | ✅ Yes   |
| `/refactor-plugin`               | Yes            | Complete plugin refactoring workflow                                                                        | ✅ Yes   |
| `/refactor-skill`                | Yes            | Split oversized skills into smaller focused skills                                                          | ✅ Yes   |
| `/start-refactor-task`           | Yes            | Execute individual refactoring tasks                                                                        | ✅ Yes   |

### Agents (6)

| Agent                      | Model  | Tools                                                                | Purpose                                                                  | Verified |
| -------------------------- | ------ | -------------------------------------------------------------------- | ------------------------------------------------------------------------ | -------- |
| `refactor-planner`         | sonnet | Read, Grep, Glob, Write                                              | Analyze plugins and create refactoring plans                             | ✅ Yes   |
| `refactor-executor`        | sonnet | Read, Write, Edit, Grep, Glob, Bash, Task                            | Execute refactoring tasks from plans                                     | ✅ Yes   |
| `refactor-validator`       | sonnet | Read, Grep, Glob, Bash                                               | Validate refactoring completeness and quality                            | ✅ Yes   |
| `subagent-refactorer`      | sonnet | Read, Write, Edit, Grep, Glob, WebFetch, WebSearch, Skill, MCP tools | Refactor Claude agents using Anthropic prompt engineering best practices | ✅ Yes   |
| `contextual-ai-documentation-optimizer` | sonnet | (inherits)                                                           | Optimize prompts, SKILL.md, CLAUDE.md for Claude comprehension           | ✅ Yes   |
| `plugin-assessor`          | sonnet | (inherits)                                                           | Analyze plugins for structure, frontmatter, and quality                  | ✅ Yes   |

### Commands (0)

### Scripts (5)

| Script                        | Purpose                                            | Verified Working                              |
| ----------------------------- | -------------------------------------------------- | --------------------------------------------- |
| `create_plugin.py`            | Interactive plugin scaffolding                     | ✅ Yes - creates .claude-plugin/, plugin.json |
| `plugin_validator.py`         | Comprehensive plugin validation with token metrics| ✅ Yes - validates plugins, skills, agents, commands |
| `validate-task-file.sh`       | Validate refactoring task file format              | ✅ Yes - bash script                          |
| `fix_tool_formats.py`         | Fix invalid tool format patterns in frontmatter    | ✅ Yes - scans ~/.claude and repos            |
| `README.md`                   | Script documentation                               | ✅ Yes - documents fix_tool_formats.py        |

---

## Consolidated Agents (v2.6.0)

**As of version 2.6.0**, the following agents were consolidated from external locations into plugin-creator to make the plugin self-contained:

### subagent-refactorer

**Source:** Previously at `~/.claude/agents/subagent-refactorer.md`
**Purpose:** Refactor Claude Code subagents using Anthropic prompt engineering best practices
**Usage:** AGENT_OPTIMIZE task types in refactoring workflows
**Key Features:**

- Researches current Anthropic documentation before refactoring
- Applies strategic XML tagging and Constitutional AI patterns
- Optimizes for Sonnet 4.5 vs Opus 4.1 model selection
- Generates analysis reports with citations and validation checklists

### contextual-ai-documentation-optimizer

**Source:** Previously at `./.claude/agents/contextual-ai-documentation-optimizer.md`
**Purpose:** Optimize prompts, SKILL.md, and CLAUDE.md files for Claude comprehension
**Usage:** DOC_IMPROVE and ORPHAN_RESOLVE task types
**Key Features:**

- Enables `prompt-optimization-claude-45` skill (external dependency)
- Applies positive framing over negative constraints
- Uses concrete examples over abstract descriptions
- Front-loads critical instructions

**Note:** Requires `prompt-optimization-claude-45` skill from separate plugin.

### plugin-assessor

**Source:** Previously at `./.claude/agents/plugin-assessor.md`
**Purpose:** Analyze plugins for structure, frontmatter, schema compliance, and quality
**Usage:** Validation tasks and pre-marketplace review
**Key Features:**

- Comprehensive reference file audit (orphan detection)
- Cross-reference validation and link graph analysis
- Frontmatter schema validation against official specs
- Generates detailed assessment reports with scoring

**Skills loaded:** `claude-skills-overview-2026`, `claude-plugins-reference-2026`, `claude-hooks-reference-2026` (all included in plugin-creator)

---

## Verified Workflows

### Workflow 1: Create New Plugin

**Entry Point:** Run script directly
**Verified:** ✅ Yes

```bash
# Interactive plugin scaffolding
uv run plugins/plugin-creator/scripts/create_plugin.py
```

**What It Does:**

1. Prompts for plugin name, description, author
2. Creates `.claude-plugin/` directory
3. Creates `plugin.json` with validated schema
4. Optionally creates `skills/`, `agents/` directories
5. Self-validates before reporting success

**Validation:** Runs `claude plugin validate` internally

---

### Workflow 2: Create New Agent

**Entry Point:** `/agent-creator` skill
**Verified:** ✅ Yes

**Trigger Phrases:**

- "Create a new agent"
- "Add an agent to {plugin}"
- "I need an agent for {task}"

**Agent Creation Process:**

1. **Discovery Phase:**

   - Reads existing agents in `.claude/agents/`
   - Identifies similar agents as templates
   - Reviews archetype templates

2. **Requirements Gathering:**

   - Uses AskUserQuestion for: purpose, triggers, tools, model, skills

3. **Template Selection:**

   - Presents options: existing project agents, role archetypes, from scratch
   - User selects via AskUserQuestion

4. **Agent File Creation:**

   - Creates frontmatter with validated fields
   - Writes agent body with workflow, quality standards

5. **Scope Determination** (uses AskUserQuestion):

   - **Project-level:** Saves to `.claude/agents/{name}.md`
   - **User-level:** Saves to `~/.claude/agents/{name}.md`
   - **Plugin:** Saves to `{plugin}/agents/{name}.md` + updates plugin.json

6. **Validation:**
   - Runs `plugin_validator.py` on agent file
   - If plugin agent: runs `claude plugin validate {plugin-path}`

**Outputs:**

- Agent file at appropriate location
- Updated plugin.json (if plugin agent)
- Validation report

---

### Workflow 3: Create New Skill

**Entry Point:** ❌ NO DEDICATED TOOL
**Verified:** ❌ GAP IDENTIFIED

**Current State:**

- No skill-creator skill exists
- No create_skill.py script exists
- Users must manually create SKILL.md files

**Gap Analysis:**

- Creating skills requires understanding frontmatter schema
- Need template selection (similar to agent-creator)
- Need scope determination (project/user/plugin)
- Need plugin.json update if plugin skill
- Need validation after creation

**Recommendation:** Create `/skill-creator` skill following agent-creator pattern

---

### Workflow 4: Create New Command

**Entry Point:** ❌ NO DEDICATED TOOL
**Verified:** ❌ GAP IDENTIFIED

**Current State:**

- No command-creator skill exists
- No create_command.py script exists
- Commands in plugins are deprecated per official docs
- Users create commands in `~/.claude/commands/` manually

**Gap Analysis:**

- Creating user-level commands requires frontmatter knowledge
- No template or validation guidance
- Commands in plugins deprecated - should document this

**Recommendation:** Create `/command-creator` skill for user-level commands only

---

### Workflow 5: Validate Plugin Components

**Entry Point:** `plugin_validator.py` script
**Verified:** ✅ Yes

**Supported File Types:**

- Complete plugins (validates all components)
- Individual SKILL.md files
- Individual agent .md files
- Individual command .md files

**Usage:**

```bash
# Validate single file
uv run plugins/plugin-creator/scripts/plugin_validator.py {path}

# Validate entire plugin directory
uv run plugins/plugin-creator/scripts/plugin_validator.py plugins/my-plugin

# Auto-fix issues
uv run plugins/plugin-creator/scripts/plugin_validator.py --fix {path}

# Validate only (no fixes)
uv run plugins/plugin-creator/scripts/plugin_validator.py --check {path}

# Verbose output with details
uv run plugins/plugin-creator/scripts/plugin_validator.py --verbose {path}

# CI mode (no color)
uv run plugins/plugin-creator/scripts/plugin_validator.py --no-color {path}
```

**What It Validates:**

- **Frontmatter schema:**
  - YAML syntax validity
  - No forbidden multiline indicators (`>-`, `|-`)
  - Required fields present
  - Field types match schema
  - Tools/skills are comma-separated strings (not YAML arrays)
- **Plugin structure:**
  - plugin.json schema compliance
  - Component path references
  - Version consistency
- **Skill complexity:**
  - Token-based complexity measurement (not line count)
  - Thresholds defined as `TOKEN_WARNING_THRESHOLD` and `TOKEN_ERROR_THRESHOLD` in `plugin_validator.py`
- **Internal links:**
  - Markdown link validity
  - Progressive disclosure structure
- **Component completeness:**
  - Required files present
  - Cross-reference validation

**What It Auto-Fixes:**

- YAML arrays → comma-separated strings
- Multiline descriptions → single-line strings
- Unquoted colons in descriptions — adds quotes to prevent YAML parsing failures
- Missing `name:` fields in plugin skills (auto-removes due to Claude Code bug)

**Schema Coverage:**

| File Type | Required Fields                   | Optional Fields                                                     | Verified |
| --------- | --------------------------------- | ------------------------------------------------------------------- | -------- |
| Skills    | None (name, description optional) | name, description, model, tools, skills, hooks, etc.                | ✅ Yes   |
| Agents    | name, description                 | model, tools, disallowedTools, permissionMode, skills, hooks, color | ✅ Yes   |
| Commands  | description                       | argument-hint, allowed-tools, model, context, agent, hooks          | ✅ Yes   |

**Error Codes:** See [ERROR_CODES.md](./scripts/ERROR_CODES.md) for complete error code reference (23 codes across 9 validators)

**SOURCE:** `plugin_validator.py` script with comprehensive validation and token-based complexity measurement

---

### Workflow 7: Validate Complete Plugin

**Entry Point:** `claude plugin validate` command
**Verified:** ✅ Yes (built-in Claude Code command)

**Usage:**

```bash
claude plugin validate {plugin-directory}
```

**What It Validates:**

- plugin.json exists in `.claude-plugin/`
- JSON syntax valid
- Required field `name` present
- `name` is kebab-case
- All paths start with `./`
- `agents` field is array of individual file paths (not directory string)
- Referenced files exist

**Common Errors:**

| Error                   | Cause                               | Fix                                              |
| ----------------------- | ----------------------------------- | ------------------------------------------------ |
| `agents: Invalid input` | Used `"./agents/"` instead of array | Change to `["./agents/file.md"]`                 |
| `name: Required`        | Missing name field                  | Add `"name": "plugin-name"`                      |
| Invalid JSON syntax     | Malformed JSON                      | Validate with `python3 -m json.tool plugin.json` |

**SOURCE:** Verified via claude-plugins-reference-2026 skill

---

### Workflow 8: Refactor Plugin

**Entry Point:** `/refactor-plugin` skill
**Verified:** ✅ Yes

**Trigger Phrases:**

- "Refactor {plugin-name} plugin"
- "Analyze plugin for refactoring"
- "Create refactoring plan for {plugin}"

**Process Flow:**

1. **Assessment** (delegates to `@"plugin-creator:refactor-planner (agent)"`):

   - Analyzes plugin structure
   - Identifies oversized skills (exceeding validator token threshold SK006/SK007)
   - Checks agent descriptions for weak triggers
   - Detects orphaned files
   - Creates assessment report

2. **Design:**

   - Creates `refactor-design-{slug}.md` with strategy

3. **Task Planning:**

   - Creates `tasks-refactor-{slug}.md` with executable tasks
   - Maps dependencies
   - Identifies parallelization opportunities

4. **Execution** (delegates to `@"plugin-creator:refactor-executor (agent)"`):

   - Executes tasks in dependency order
   - Runs parallel where possible
   - Tracks completion status

5. **Validation** (delegates to `@"plugin-creator:refactor-validator (agent)"`):
   - Verifies refactoring achieved goals
   - Checks for regressions
   - Creates follow-up tasks if issues found

**Task Types Handled:**

| Type              | Handler                                 | Verified |
| ----------------- | --------------------------------------- | -------- |
| `SKILL_SPLIT`     | `/plugin-creator:refactor-skill` skill  | ✅ Yes   |
| `AGENT_OPTIMIZE`  | `subagent-refactorer` agent             | ✅ Yes   |
| `DOC_IMPROVE`     | `contextual-ai-documentation-optimizer` agent        | ✅ Yes   |
| `DOC_UPDATER_ADD` | `/plugin-creator:add-doc-updater` skill | ✅ Yes   |
| `ORPHAN_RESOLVE`  | Manual or context optimizer             | ✅ Yes   |
| `STRUCTURE_FIX`   | Direct implementation                   | ✅ Yes   |

**When DOC_UPDATER_ADD is identified:**
- Skill wraps external documentation (API specs, CLI refs, frameworks)
- Documentation source updates regularly (weekly, monthly)
- Skill would benefit from automated sync vs. manual updates

**OUTPUT:** Task files in `.claude/plan/` directory

---

### Workflow 9: Split Oversized Skill

**Entry Point:** `/refactor-skill` skill
**Verified:** ✅ Yes

**Model:** opus (requires complex reasoning)

**Trigger Phrases:**

- "Split the {skill-name} skill"
- "Refactor oversized skill"
- "This skill exceeds the token limit"

**Process:**

1. Reads existing skill SKILL.md
2. Identifies logical boundaries and domains
3. Designs split plan
4. Generates new SKILL.md files for each extracted skill
5. Validates 100% content migration (no loss)
6. Creates cross-references between new skills
7. Updates original skill as facade/meta-skill

**Quality Requirements:**

- No content loss
- Complete fidelity
- Backwards compatibility maintained
- All cross-references valid

---

### Workflow 12: Add Documentation Updater to Skill

**Entry Point:** `/add-doc-updater` skill
**Verified:** ✅ Yes

**Trigger Phrases:**

- "Add doc sync to {skill}"
- "Automate documentation updates for {skill}"
- "This skill needs to wrap {external docs}"

**Process Flow:**

1. **Phase 0 - Variable Collection:**

   - Validate target skill path
   - Infer defaults (skill name, cooldown 7 days)
   - Collect 6 template variables via AskUserQuestion
   - Confirm before proceeding

2. **Phase 1 - Implementation:**

   - Load `/python3-development` orchestration
   - Substitute template variables
   - Delegate to `@python-cli-architect` agent
   - Creates `scripts/update-{LOCAL_DOC_DIR}-docs.py`

3. **Phase 2 - Code Review:**

   - Delegate to `@python-code-reviewer` agent
   - Validate: ReDoS-safe regex, atomic operations, link transformation
   - Loop back to Phase 1 on failure

4. **Phase 3 - Quality Gates:**

   - Sequential: ruff format → ruff check → mypy → pyright → prek
   - All must pass to proceed
   - Loop back to Phase 1 on failure

5. **Phase 4 - Testing & Validation:**

   - 7-point checklist: execution, file existence, Hugo shortcode removal, link sampling, SKILL.md integration, cooldown enforcement, force flag
   - Loop back to Phase 1 on failure

6. **Phase 5 - Integration:**
   - Update SKILL.md Execution Protocol section
   - Add .gitignore entries for lock files and downloaded docs
   - Integration test via general-purpose agent
   - Loop back to Phase 1 if agent cannot discover documentation

**Success Criteria:**

- Script passes all quality gates
- All 7 validation tests pass
- SKILL.md includes execution protocol
- Integration test demonstrates autonomous documentation access
- Cooldown prevents excessive upstream requests
- Force flag allows manual override

**OUTPUT:** Documentation updater ready for use

**SOURCE:** Verified by reading add-doc-updater/SKILL.md lines 1-354

---

### Workflow 11: Fix Tool Formatting Issues

**Entry Point:** `fix_tool_formats.py` script
**Verified:** ✅ Yes

**Purpose:** Fix invalid tool format patterns in frontmatter across entire codebase

**Usage:**

```bash
# Scans ~/.claude and ~/repos recursively
uv run plugins/plugin-creator/scripts/fix_tool_formats.py
```

**What It Fixes:**

1. **YAML list → comma-separated string:**

   ```yaml
   # Before
   tools:
     - Read
     - Grep

   # After
   tools: Read, Grep
   ```

2. **JSON array → comma-separated string:**

   ```yaml
   # Before
   tools: ["Read", "Grep"]

   # After
   tools: Read, Grep
   ```

**Scan Locations:**

- `~/.claude/agents/**/*.md`
- `~/.claude/commands/**/*.md`
- `~/.claude/skills/**/SKILL.md`
- `~/repos/**/.claude/**` (all markdown in .claude directories)

**Why This Matters:** Invalid formats become "evidence" in future Grep searches, creating feedback loop where AI learns incorrect patterns from its own mistakes.

**SOURCE:** Verified by reading `scripts/README.md` and `fix_tool_formats.py`

---

## Identified Gaps

### Gap 1: No Skill Creator

**Current State:** Can create agents, cannot create skills
**Impact:** Users must manually create SKILL.md files without guidance
**Recommendation:** Create `/skill-creator` skill similar to `/agent-creator`

**Required Capabilities:**

- Template selection (general-purpose, domain-specific, etc.)
- Frontmatter generation with validation
- Scope determination (project/user/plugin)
- Plugin.json update if plugin skill
- Post-creation validation

---

### Note: Skills Are Commands

**Observation from official docs (line 9 of claude-skills-overview-2026):**

> "Skills and slash commands are now unified - they are the same system. A file at `.claude/commands/review.md` and a skill at `.claude/skills/review/SKILL.md` both create `/review` and work identically. Skills are the recommended approach."

**Conclusion:** No separate command-creator needed. Creating a skill with `user-invocable: true` creates a slash command.

---

### Gap 2: Incomplete Validation Coverage

**Validated:**

- ✅ Frontmatter schema (via plugin_validator.py)
- ✅ Plugin.json syntax (via claude plugin validate)
- ✅ Skill complexity (token-based, via plugin_validator.py)

**Not Validated:**

- ❌ Agent frontmatter vs skill frontmatter differences (plugin_validator.py validates both but does not enforce agent-specific constraints)
- ❌ Command frontmatter (no dedicated validator beyond general plugin_validator.py)
- ❌ Cross-references between components (e.g., agent references non-existent skill)

**Recommendation:** Extend validation to cover agent-specific checks and cross-reference validation

---

## Usage Patterns for the Model

### Pattern 1: User Wants to Create Something

**If user says:** "Create a new {plugin|agent|skill|command}"

```
1. Identify what they want to create
2. Check if tool exists:
   - Plugin → uv run scripts/create_plugin.py
   - Agent → /agent-creator
   - Skill → ❌ GAP - no tool (manual creation)
   - Command → ❌ GAP - no tool (manual creation)
3. If gap: Offer manual creation with validation guidance
```

---

### Pattern 2: User Wants to Validate

**If user says:** "Validate {frontmatter|plugin|skill}"

```
1. Determine validation type:
   - Frontmatter → uv run scripts/plugin_validator.py validate {path}
   - Complete plugin → uv run scripts/plugin_validator.py {plugin-path} then claude plugin validate {path}
   - Skill structure → uv run scripts/plugin_validator.py {path}
2. Run appropriate validator
3. Report results with specific file:line references
```

---

### Pattern 3: User Wants to Refactor

**If user says:** "Refactor {plugin}" or "This skill is too large"

```
1. If skill may be oversized:
   - Run uv run plugins/plugin-creator/scripts/plugin_validator.py {path} to check token complexity
   - If validator reports SK006 or SK007: Offer /refactor-skill for individual skill or /refactor-plugin for whole plugin
2. If general plugin refactoring:
   - Delegate to refactor-planner agent
   - Review generated plan with user
   - Delegate execution to refactor-executor agent
   - Validate with refactor-validator agent
```

---

### Pattern 4: User Reports Validation Errors

**If user says:** "Getting 'agents: Invalid input'" or similar

```
1. Read their plugin.json
2. Check against common errors:
   - agents field must be array of file paths
   - tools fields must be comma-separated strings
   - All paths must start with ./
3. Offer fix_tool_formats.py for tool format issues
4. Offer plugin_validator.py fix for frontmatter issues
```

---

## Script Execution Paths

All scripts use absolute paths or are executable with `uv run`:

```bash
# From anywhere in the repository:
uv run plugins/plugin-creator/scripts/plugin_validator.py validate {path}
uv run plugins/plugin-creator/scripts/create_plugin.py

# Scripts that work from their directory:
cd plugins/plugin-creator/scripts
./validate-task-file.sh {path}
```

**Important:** Scripts expect to be run from repository root or use `${CLAUDE_PLUGIN_ROOT}` variable.

---

## Plugin System Fundamentals

### Plugin Caching and File Resolution

**CRITICAL:** Claude Code copies plugins to a cache directory rather than using them in-place.

**How it works:**

- Marketplace plugins with relative paths: The `source` path is copied recursively
- Plugins with `.claude-plugin/plugin.json`: The directory containing `.claude-plugin/` is copied recursively

**Path traversal limitations:**

- Plugins CANNOT reference files outside their directory (`../shared-utils` will FAIL after installation)
- External files are NOT copied to the cache

**Solutions for external dependencies:**

1. **Use symlinks:** Create symlinks within plugin directory (symlinks are followed during copy)

   ```bash
   ln -s /path/to/shared-utils ./shared-utils
   ```

2. **Restructure marketplace:** Set source to parent directory that contains all required files

**SOURCE:** Lines 350-398 of claude-plugins-reference-2026/SKILL.md

### Installation Scopes

When installing plugins, the scope determines where the plugin is available:

| Scope     | Settings file                 | Use case                                                 |
| --------- | ----------------------------- | -------------------------------------------------------- |
| `user`    | `~/.claude/settings.json`     | Personal plugins available across all projects (default) |
| `project` | `.claude/settings.json`       | Team plugins shared via version control                  |
| `local`   | `.claude/settings.local.json` | Project-specific plugins, gitignored                     |
| `managed` | `managed-settings.json`       | Managed plugins (read-only, update only)                 |

**Usage examples:**

```bash
# Install to user scope (default)
claude plugin install formatter@my-marketplace

# Install to project scope (shared with team)
claude plugin install formatter@my-marketplace --scope project

# Install to local scope (gitignored)
claude plugin install formatter@my-marketplace --scope local
```

**SOURCE:** Lines 401-410 of claude-plugins-reference-2026/SKILL.md

### Environment Variables

When commands execute, they have access to:

| Variable                | Value                                                  | Usage                                 |
| ----------------------- | ------------------------------------------------------ | ------------------------------------- |
| `${CLAUDE_PLUGIN_ROOT}` | Absolute path to plugin directory                      | Used in commands to reference scripts |
| `${CLAUDE_PROJECT_DIR}` | Project root directory (where Claude Code was started) | Project-relative paths                |
| `$ARGUMENTS`            | Command arguments from user                            | Passed to command's bash execution    |

**SOURCE:** Lines 415-434 of claude-plugins-reference-2026/SKILL.md

---

## Quality Standards Enforced

### Skill Size Limits

Run `uv run plugins/plugin-creator/scripts/plugin_validator.py <skill-path>` after writing and follow its guidance. Thresholds are defined as `TOKEN_WARNING_THRESHOLD` (SK006) and `TOKEN_ERROR_THRESHOLD` (SK007) in `plugin_validator.py` — not as line counts. SK006 triggers `references/` extraction; SK007 requires skill splitting.

### Frontmatter Requirements

**Skills:**

- `name`: Optional (uses directory name if omitted)
- `description`: Optional (uses first paragraph if omitted)

**Agents:**

- `name`: Required (lowercase, hyphens, max 64 chars)
- `description`: Required (trigger keywords, max 1024 chars)
- `model`: Must be sonnet/opus/haiku/inherit if specified
- `tools`: Must be comma-separated string (not YAML array)

**Commands:**

- `description`: Required
- `allowed-tools`: Must be comma-separated string (not YAML array)

### Plugin.json Requirements

**Required fields:**

- `name`: Required, kebab-case

**Component path fields:**

| Field          | Type           | Description                                         | Example                                  |
| -------------- | -------------- | --------------------------------------------------- | ---------------------------------------- |
| `commands`     | string\|array  | Additional command files/directories                | `"./custom/cmd.md"` or `["./cmd1.md"]`   |
| `agents`       | string\|array  | Additional agent files or directories               | `"./custom/agents/"` or `["./agent.md"]` |
| `skills`       | string\|array  | Additional skill directories                        | `"./custom/skills/"`                     |
| `hooks`        | string\|object | Hook config path or inline config                   | `"./hooks.json"`                         |
| `mcpServers`   | string\|object | MCP config path or inline config                    | `"./mcp-config.json"`                    |
| `outputStyles` | string\|array  | Additional output style files/directories           | `"./styles/"`                            |
| `lspServers`   | string\|object | Language Server Protocol config (code intelligence) | `"./.lsp.json"`                          |

**Path behavior rules:**

- Custom paths supplement default directories - they don't replace them
- If `commands/` exists, it's loaded in addition to custom command paths
- All paths must be relative and start with `./`
- Multiple paths can be specified as arrays

**SOURCE:** Lines 75-91 of claude-plugins-reference-2026/SKILL.md

---

## Verification Status

This documentation was created 2026-01-28 by:

1. Reading all skill SKILL.md files
2. Reading all agent frontmatter
3. Reading all command files
4. Reading all script files
5. Testing scripts where possible
6. Cross-referencing with official Claude Code documentation

**Verification Method:** Direct file reading and execution verification, not assumption-based.

**Gaps Explicitly Identified:** Listed in "Identified Gaps" section above.

---

## Script Consolidation Recommendation

**Current State:** Validation and linting functionality is spread across multiple scripts:

- `plugin_validator.py` - Frontmatter schema validation, skill complexity (token-based), internal links
- `validate-task-file.sh` - Task file validation (bash)
- `fix_tool_formats.py` - Tool format fixing

**Recommendation:** Consolidate into a single cross-platform Python script: `lint-claude-plugin.py`

**Requirements:**

1. **Single Script:** Combine all validation functionality into one Python 3.11+ script
2. **Cross-Platform:** Pure Python, no bash dependencies, works on Windows/Linux/macOS
3. **Pre-commit Integration:** Compatible with `.pre-commit-config.yaml` hooks
4. **Token-Based Metrics:** Use `tiktoken` library to measure skill complexity in tokens, not lines
   - Line count is a poor proxy for complexity
   - Token count directly measures what Claude processes
   - Thresholds should be token-based (e.g., warn at X tokens, error at Y tokens)
5. **Unified Validation:** Single entry point that validates:
   - Frontmatter schema (skills, agents, commands)
   - Plugin.json structure
   - Skill complexity via token count
   - Internal link validity
   - Progressive disclosure structure
   - Tool format correctness

**Benefits:**

- Single tool to install and maintain
- Consistent behavior across platforms
- Better complexity measurement (tokens vs lines)
- Pre-commit hook compatibility
- Reduced maintenance burden

**Implementation Priority:** Medium - Current scripts work but consolidation would improve usability

---

## LSP Server Integration

Plugins can provide Language Server Protocol (LSP) servers for real-time code intelligence:

**Capabilities:**

- **Instant diagnostics:** Claude sees errors and warnings immediately after each edit
- **Code navigation:** go to definition, find references, hover information
- **Language awareness:** type information and documentation for code symbols

**Configuration format (`.lsp.json` or inline in `plugin.json`):**

```json
{
  "lspServers": {
    "python": {
      "command": "pyright-langserver",
      "args": ["--stdio"],
      "extensionToLanguage": {
        ".py": "python"
      }
    }
  }
}
```

**Required fields:**

| Field                 | Description                                  |
| --------------------- | -------------------------------------------- |
| `command`             | The LSP binary to execute (must be in PATH)  |
| `extensionToLanguage` | Maps file extensions to language identifiers |

**Optional fields:**

| Field                   | Description                                               |
| ----------------------- | --------------------------------------------------------- |
| `args`                  | Command-line arguments for the LSP server                 |
| `transport`             | Communication transport: `stdio` (default) or `socket`    |
| `env`                   | Environment variables to set when starting the server     |
| `initializationOptions` | Options passed to the server during initialization        |
| `settings`              | Settings passed via `workspace/didChangeConfiguration`    |
| `workspaceFolder`       | Workspace folder path for the server                      |
| `startupTimeout`        | Max time to wait for server startup (milliseconds)        |
| `shutdownTimeout`       | Max time to wait for graceful shutdown (milliseconds)     |
| `restartOnCrash`        | Whether to automatically restart the server if it crashes |
| `maxRestarts`           | Maximum number of restart attempts before giving up       |

**CRITICAL:** LSP servers require separate binary installation. LSP plugins configure Claude Code's connection to a language server but don't include the server itself.

**Available LSP plugins:**

| Plugin           | Language server  | Install command                                                                            |
| ---------------- | ---------------- | ------------------------------------------------------------------------------------------ |
| `pyright-lsp`    | Pyright (Python) | `pip install pyright` or `npm install -g pyright`                                          |
| `typescript-lsp` | TypeScript LS    | `npm install -g typescript-language-server typescript`                                     |
| `rust-lsp`       | rust-analyzer    | See [rust-analyzer installation](https://rust-analyzer.github.io/manual.html#installation) |

**SOURCE:** Lines 271-347 of claude-plugins-reference-2026/SKILL.md

---

## CLI Commands Reference

The model MUST use these CLI commands for plugin management:

**Install plugin:**

```bash
claude plugin install <plugin> [--scope user|project|local]
```

**Uninstall plugin:**

```bash
claude plugin uninstall <plugin> [--scope user|project|local]
# Aliases: remove, rm
```

**Enable/disable plugin:**

```bash
claude plugin enable <plugin> [--scope user|project|local]
claude plugin disable <plugin> [--scope user|project|local]
```

**Update plugin:**

```bash
claude plugin update <plugin> [--scope user|project|local|managed]
```

**Validate plugin:**

```bash
claude plugin validate <plugin-directory>
/plugin validate <plugin-directory>  # In Claude Code session
```

**Testing without installation (session only):**

```bash
# Load plugin for current session only
claude --plugin-dir ./my-plugin

# Load multiple plugins
claude --plugin-dir ./plugin-one --plugin-dir ./plugin-two
```

**SOURCE:** Lines 565-729 of claude-plugins-reference-2026/SKILL.md

---

## Next Desired Outcomes (Pre-Context Gathering Phase)

This section lists discrete desired outcomes ready for Stage 2 (Planning with RT-ICA) following the Stateless Software Engineering Framework. Each outcome can be worked through independently by specialist agents.

**Methodology Reference**: [stateless-software-engineering-framework.md](https://github.com/bitflight-devops/stateless-agent-methodology/blob/main/stateless-software-engineering-framework.md)

**SSE Stage Progression for Each Outcome**:

1. **Stage 1: Discovery** (COMPLETE - documented below)
2. **Stage 2: Planning** with RT-ICA (NEXT - requires specialist agent)
3. **Stage 3: Context Integration** (validates plan against codebase)
4. **Stage 4: Task Decomposition** (creates executable tasks)
5. **Stage 5: Execution** (stateless implementation)
6. **Stage 6: Forensic Review** (validates completion)
7. **Stage 7: Final Verification** (certifies outcome achievement)

---

### Outcome 1: User-Invocable `/create-plugin` Workflow

**Objective**: Provide systematic plugin creation workflow following SSE methodology

**Current State** (Stage 1 Discovery Output):

- Script exists: `scripts/create_plugin.py` (interactive CLI)
- Skill exists: `plugin-creator` (delegates to agents)
- Gap: No RT-ICA phase, no quality gates, script not integrated with agent orchestration

**Prerequisites to Verify (RT-ICA Inputs)**:

- Existing agent orchestration patterns in plugin-creator
- Location of RT-ICA skill/methodology
- Template structure for user-invocable skills
- Validation scripts that must run at each gate
- Integration points with existing create_plugin.py script

**Dependencies on Existing Components**:

- `claude-plugins-reference-2026` skill (reference documentation)
- `claude-skills-overview-2026` skill (skill structure guidance)
- `plugin-creator` orchestrator (existing delegation pattern)
- `plugin_validator.py` script (quality gate)
- `claude plugin validate` command (structure validation)

**Success Criteria**:

- Skill file created at `skills/create-plugin/SKILL.md`
- Workflow follows 6-phase structure from existing `plugin-dev:create-plugin` skill
- RT-ICA phase explicitly included before planning
- Quality gates run after each phase
- Backwards compatible with existing create_plugin.py script
- Validates using existing validation tools

**Acceptance Criteria**:

- `/create-plugin` command appears in skill list
- Triggering skill starts Discovery phase with RT-ICA
- Each phase blocks on validation gate before proceeding
- Final output is validated plugin directory structure
- Documentation updated in CLAUDE.md

**Known Unknowns** (requires RT-ICA assessment):

- How to integrate script-based scaffolding with skill-based orchestration?
- Should RT-ICA be embedded in skill or delegated to separate skill?
- What validation gates block progression between phases?
- How to handle user approval checkpoints?

---

### Outcome 2: User-Invocable `/extend-plugin` Workflow

**Objective**: Add new components to existing plugins with validation and integration checks

**Current State** (Stage 1 Discovery Output):

- No dedicated workflow exists
- Users manually create components without guidance
- No validation of compatibility with existing structure
- No automatic plugin.json updates

**Prerequisites to Verify (RT-ICA Inputs)**:

- Existing plugin.json update patterns
- Component type detection logic
- Validation requirements for each component type
- Integration test patterns for new components
- Scope determination patterns (project/user/plugin)

**Dependencies on Existing Components**:

- `skill-creator` skill (existing)
- `agent-creator` skill (existing)
- `write-frontmatter-description` skill (description generation)
- `plugin_validator.py` script (validation gate)
- `claude plugin validate` command (structure validation)

**Success Criteria**:

- Skill file created at `skills/extend-plugin/SKILL.md`
- Workflow analyzes existing plugin before adding components
- Follows existing naming and structure conventions
- Updates plugin.json automatically
- Validates new component against plugin schema

**Acceptance Criteria**:

- `/extend-plugin <plugin-path> --add <type>` command works
- Detects and follows existing plugin patterns
- Updates plugin.json with new component reference
- Runs validation suite on modified plugin
- Reports integration status (COMPATIBLE / CONFLICTS)

**Known Unknowns** (requires RT-ICA assessment):

- How to detect existing naming conventions automatically?
- What constitutes a "conflict" with existing structure?
- Should this skill create missing skill-creator if needed?
- How to handle plugins with custom directory structures?

---

### Outcome 3: User-Invocable `/validate-plugin` Workflow

**Objective**: Comprehensive validation suite with quality scoring and reporting

**Current State** (Stage 1 Discovery Output):

- Partial validation: `plugin_validator.py` (frontmatter, complexity, links)
- Partial validation: `claude plugin validate` (structure only)
- No unified entry point, no quality scoring, no marketplace readiness check

**Prerequisites to Verify (RT-ICA Inputs)**:

- All validation scripts and their exit codes
- Marketplace submission requirements
- Quality scoring rubric (if exists)
- Report generation format
- Integration with existing plugin-assessor agent

**Dependencies on Existing Components**:

- `plugin_validator.py` script (schema validation, complexity, links)
- `claude plugin validate` CLI (structure validation)
- `plugin-assessor` agent (comprehensive analysis)
- Reference skills for validation rules

**Success Criteria**:

- Skill file created at `skills/validate-plugin/SKILL.md`
- Runs all existing validators in sequence
- Generates unified report with pass/fail per category
- Includes marketplace readiness score
- Validates version consistency across files

**Acceptance Criteria**:

- `/validate-plugin <plugin-path>` command works
- Creates `validation-report.md` artifact
- Reports scores: structure, frontmatter, quality, marketplace
- Exit code 0 = pass, non-zero = fail
- Can run in CI/pre-commit context

**Known Unknowns** (requires RT-ICA assessment):

- What is the complete marketplace submission checklist?
- How to score quality (weighted vs binary)?
- Should validation auto-fix issues or only report?
- Integration with existing plugin-assessor agent workflow?

---

### Outcome 4: ~~Create `/skill-creator` Skill~~ (ALREADY EXISTS - REMOVED)

**Status**: This outcome has been removed from the roadmap.

**Correction**: During Stage 1 Discovery, the skill-creator skill was incorrectly identified as a gap. The skill already exists at `plugins/plugin-creator/skills/skill-creator/SKILL.md` and is fully functional.

**Verification**: See Component Inventory line 46 in this document which lists `/skill-creator` as an existing user-invocable skill.

**Impact on Outcome 2**: Outcome 2 (/extend-plugin workflow) listed skill-creator as a dependency. That dependency is satisfied - no blocking issue exists.

---

### Outcome 5: Pre-Commit Validation Hook

**Objective**: Automatically validate plugins on git commit to prevent invalid state

**Current State** (Stage 1 Discovery Output):

- No pre-commit hook exists for plugin validation
- Repository uses `prek` (Rust pre-commit replacement)
- Validation scripts exist but aren't hooked into git workflow

**Prerequisites to Verify (RT-ICA Inputs)**:

- Existing `.pre-commit-config.yaml` structure
- `prek` hook configuration format
- Changed file detection logic
- Performance constraints (hook must be fast)
- Error reporting format for git hooks

**Dependencies on Existing Components**:

- `.pre-commit-config.yaml` (configuration file)
- `prek` tool (pre-commit runner)
- `plugin_validator.py` script (validation logic, quality checks)
- `claude plugin validate` command (structure validation)

**Success Criteria**:

- Hook configuration added to `.pre-commit-config.yaml`
- Hook detects changed plugins from git diff
- Runs only validators relevant to changed files
- Fast execution (<5s for typical changes)
- Clear error messages with file:line references

**Acceptance Criteria**:

- Pre-commit hook blocks commits with validation errors
- Hook passes when no plugin files changed
- Hook validates only changed plugins (not entire repo)
- Error output shows which validator failed and why
- Hook can be bypassed with `--no-verify` if needed

**Known Unknowns** (requires RT-ICA assessment):

- Does `prek` use identical syntax to `pre-commit`?
- What is acceptable performance threshold?
- Should hook auto-fix or only validate?
- How to handle plugin.json changes vs component changes?
- Integration with existing markdown linting hooks?

---

### Outcome 6: Automated Version Bumping

**Objective**: Automatically increment plugin and marketplace versions when plugins change

**Current State** (Stage 1 Discovery Output):

- No automatic version bumping exists
- Plugin versions manually updated in plugin.json
- Marketplace version manually updated in marketplace.json
- Users receive stale versions when pulling updates

**Prerequisites to Verify (RT-ICA Inputs)**:

- Semantic versioning policy (major/minor/patch rules)
- Git commit message conventions for version type detection
- File patterns that trigger version bumps
- Version field locations (plugin.json, marketplace.json)
- Existing version bumping tools/libraries

**Dependencies on Existing Components**:

- `.pre-commit-config.yaml` (hook integration)
- `prek` tool (pre-commit runner)
- `plugin.json` files (version storage)
- `.claude-plugin/marketplace.json` (marketplace version)
- Git commit messages (version type signals)

**Success Criteria**:

- Hook configuration added to `.pre-commit-config.yaml`
- Detects changed plugins and bumps their versions
- Bumps marketplace version when any plugin changes
- Follows semantic versioning rules
- Stages updated version files automatically

**Acceptance Criteria**:

- Plugin version incremented when plugin files change
- Marketplace version incremented when plugins change
- Version bump type determined from commit message or change type
- Updated JSON files staged for commit automatically
- Version history maintained in CLAUDE.md

**Status**: ✅ **IMPLEMENTED** (2026-01-29)

**Implementation**:

- Script: `plugins/plugin-creator/scripts/auto_sync_manifests.py`
- Hook: Configured in `.pre-commit-config.yaml` as `auto-sync-manifests`
- Documentation: `plugins/plugin-creator/scripts/README-auto-sync.md`

**How It Works**:

The pre-commit hook automatically:

1. Detects CRUD operations on plugins and components from git staged changes
2. Updates `plugin.json` component arrays (skills, agents, commands) with proper `./` paths
3. Bumps plugin versions following semantic versioning:
   - **Major** (X.0.0): Component deleted (breaking change)
   - **Minor** (0.X.0): Component added (new feature)
   - **Patch** (0.0.X): Component modified
4. Updates `marketplace.json` plugin registry
5. Bumps marketplace version:
   - **Major** (X.0.0): Plugin deleted
   - **Minor** (0.X.0): Plugin added
   - **Patch** (0.0.X): Plugin component changed
6. Stages updated manifest files automatically

**Version Bumping Rules**:

```
Component Change → Plugin Version    → Marketplace Version
─────────────────────────────────────────────────────────
+ New skill      → Minor (0.1.0→0.2.0) → Patch (1.0.0→1.0.1)
- Delete agent   → Major (0.1.0→1.0.0) → Patch (1.0.0→1.0.1)
~ Modify command → Patch (0.1.0→0.1.1) → Patch (1.0.0→1.0.1)

+ New plugin     → N/A                 → Minor (1.0.0→1.1.0)
- Delete plugin  → N/A                 → Major (1.0.0→2.0.0)
```

**Usage**:

Automatic during `git commit`:

```bash
git add plugins/my-plugin/skills/new-skill/SKILL.md
git commit -m "feat(my-plugin): add new skill"
# ✅ Updated my-plugin → 1.1.0 (+1 components)
# ✅ Updated marketplace → 2.0.1
```

Manual execution:

```bash
uv run -q --no-sync plugins/plugin-creator/scripts/auto_sync_manifests.py
```

**Resolution of Known Unknowns**:

- ✅ Major/minor/patch determined by file operation type (delete/add/modify)
- ✅ Version type auto-detected from git diff, independent of commit message
- ✅ Multiple plugins: Each bumps independently, marketplace bumps once
- ✅ Version conflicts: Pre-commit stage prevents concurrent modification
- ✅ Schema integration: Validated against marketplace.json structure

---

### Outcome 7: Consolidated Cross-Platform Linter (`lint-claude-plugin.py`)

**Objective**: Replace bash validation scripts with unified Python linter using token-based complexity metrics

**Current State** (Stage 1 Discovery Output):

- Validation spread across multiple scripts (Python + Bash)
- Line count metrics don't reflect actual complexity
- Bash scripts not cross-platform compatible
- No token-based complexity measurement

**Prerequisites to Verify (RT-ICA Inputs)**:

- Token counting library (`tiktoken`) API and usage
- Token thresholds for skill complexity warnings
- All validation rules from existing scripts
- Pre-commit integration requirements
- Output format for CI/IDE integration

**Dependencies on Existing Components**:

- `plugin_validator.py` (frontmatter validation, complexity, quality check logic)
- `validate-task-file.sh` (task file validation logic)
- `fix_tool_formats.py` (auto-fix logic)

**Success Criteria**:

- Single Python script at `scripts/lint-claude-plugin.py`
- Validates frontmatter, structure, complexity, links, tool formats
- Uses token count instead of line count for complexity
- Works on Windows/Linux/macOS without bash
- Can auto-fix issues with `--fix` flag

**Acceptance Criteria**:

- All validation logic from bash scripts ported to Python
- Token-based complexity measurement implemented
- Exit code 0 = pass, non-zero = fail with details
- Compatible with pre-commit hooks
- Performance comparable to existing scripts (<5s typical)

**Known Unknowns** (requires RT-ICA assessment):

- What are appropriate token thresholds (warn/error)?
- Which tiktoken encoding to use (cl100k_base)?
- How to measure token count including frontmatter?
- Should auto-fix be opt-in or opt-out?
- Backwards compatibility with existing script invocations?

---

## Development Roadmap (Stage 1 Discovery - Historical)

### Missing User-Invocable Workflows

The plugin currently has workflows for refactoring but lacks streamlined user-invocable skills for core plugin lifecycle operations. Following the [SAM methodology](https://github.com/bitflight-devops/stateless-agent-methodology) processes, these workflows should be added:

#### 1. `/create-plugin` - Complete Plugin Creation Workflow

**Status:** ⚠️ Partial - has `create_plugin.py` script but no orchestrated skill

**Current State:**

- Script exists: `scripts/create_plugin.py` (interactive CLI)
- Skill exists: `plugin-creator` (delegates to agents but not systematic workflow)

**Gap:**

- No structured workflow following SAM methodology process
- No RT-ICA (Reverse Thinking - Information Completeness Assessment) phase
- No integration with validation and quality gates
- Script is standalone, not integrated with agent orchestration

**Proposed Implementation:**

```yaml
---
name: create-plugin
description: 'Guide plugin creation from concept to validated implementation following systematic methodology. Handles requirement gathering, component design, implementation, validation, and marketplace preparation. Use when starting a new plugin project.'
user-invocable: true
model: sonnet
context: fork
---

Workflow phases:
1. Discovery - RT-ICA to gather complete requirements
2. Component Planning - Determine skills/agents/commands/hooks needed
3. Structure Creation - Scaffold plugin directory and manifest
4. Implementation - Create components with validation
5. Quality Assurance - Run validators and quality checks
6. Marketplace Preparation - Prepare for distribution
```

**Dependencies:**

- Load `claude-plugins-reference-2026` skill
- Load `claude-skills-overview-2026` skill
- Delegate to `plugin-creator` orchestrator
- Run validation scripts after each phase

#### 2. `/extend-plugin` - Add Components to Existing Plugin

**Status:** ❌ Missing

**Purpose:** Add new skills, agents, commands, or hooks to an existing plugin

**Gap:**

- No systematic workflow for extending plugins
- Users manually create components without guidance
- No validation of compatibility with existing plugin structure
- No automatic plugin.json updates

**Proposed Implementation:**

```yaml
---
name: extend-plugin
description: 'Add new components (skills, agents, commands, hooks) to existing plugins with validation and integration checks. Ensures new components follow plugin conventions and update plugin.json correctly. Use when adding functionality to existing plugins.'
user-invocable: true
argument-hint: <plugin-path> --add <skill|agent|command|hook>
model: sonnet
---

Workflow phases:
1. Plugin Analysis - Read existing plugin structure and conventions
2. Component Type Selection - Ask user what to add
3. Requirement Gathering - RT-ICA for new component
4. Implementation - Create component following existing patterns
5. Integration - Update plugin.json, validate references
6. Validation - Run validators on modified plugin
```

**Dependencies:**

- Load `skill-creator`, `agent-creator` as needed
- Load `write-frontmatter-description` for descriptions
- Run `plugin_validator.py` after creation
- Run `claude plugin validate` for structure

#### 3. `/validate-plugin` - Comprehensive Plugin Validation

**Status:** ⚠️ Partial - has validation scripts but no orchestrated workflow

**Current State:**

- Script: `plugin_validator.py` (frontmatter, complexity, internal links)
- CLI: `claude plugin validate` (structure only)

**Gap:**

- No single entry point for complete validation
- No quality scoring or reporting
- No pre-commit integration
- No marketplace readiness check

**Proposed Implementation:**

```yaml
---
name: validate-plugin
description: 'Run comprehensive validation suite on plugin including structure, frontmatter, quality checks, and marketplace readiness. Generates validation report with scores and recommendations. Use before commits, releases, or marketplace submission.'
user-invocable: true
argument-hint: <plugin-path>
model: sonnet
---

Validation suite:
1. Structure validation - plugin.json schema, paths, references
2. Frontmatter validation - all SKILL.md and agent files
3. Quality checks - skill token complexity, orphaned files, broken links
4. Marketplace readiness - README, LICENSE, keywords, description
5. Version consistency - check all version fields match
6. Generate report - validation-report.md with scores and recommendations
```

**Dependencies:**

- `plugin_validator.py` (all components, skill quality checks)
- `claude plugin validate` (structure)
- `plugin-assessor` agent (comprehensive analysis)

#### 4. Pre-Commit Hook for Plugin Validation

**Status:** ❌ Missing

**Purpose:** Automatically validate plugins on git commit

**Gap:**

- No pre-commit hook for plugin validation
- Changes to plugins aren't validated before commit
- Easy to commit invalid plugin.json or frontmatter

**Proposed Implementation:**

Create `.pre-commit-config.yaml` entry:

```yaml
  - repo: local
    hooks:
      - id: validate-plugins
        name: Validate Claude Code Plugins
        entry: scripts/validate-changed-plugins.sh
        language: script
        files: '^plugins/.*/(\.claude-plugin/plugin\.json|.*\.md)$'
        pass_filenames: false
```

Create `scripts/validate-changed-plugins.sh`:

```bash
#!/usr/bin/env bash
# Validate only plugins that have changes

set -euo pipefail

# Find changed plugin directories
changed_plugins=$(git diff --cached --name-only | grep '^plugins/' | cut -d'/' -f1-2 | sort -u)

if [ -z "$changed_plugins" ]; then
  exit 0
fi

for plugin_dir in $changed_plugins; do
  echo "Validating $plugin_dir..."

  # Run frontmatter validation
  uv run plugins/plugin-creator/scripts/plugin_validator.py "$plugin_dir"

  # Run plugin structure validation
  claude plugin validate "$plugin_dir"

done

echo "✅ All changed plugins validated successfully"
```

#### 5. Automated Version Bumping

**Status:** ❌ Missing

**Purpose:** Automatically bump plugin versions and marketplace version on changes

**Gap:**

- Plugin versions not incremented when plugins change
- Marketplace version not updated automatically
- Users get stale versions when pulling updates

**Proposed Implementation:**

**Pre-commit hook for version bumping:**

```yaml
  - repo: local
    hooks:
      - id: bump-plugin-versions
        name: Bump Plugin and Marketplace Versions
        entry: scripts/bump-plugin-versions.sh
        language: script
        files: '^plugins/.*'
        pass_filenames: false
```

**Script: `scripts/bump-plugin-versions.sh`:**

```bash
#!/usr/bin/env bash
# Automatically bump versions for changed plugins

set -euo pipefail

# Find changed plugins
changed_plugins=$(git diff --cached --name-only | grep '^plugins/' | cut -d'/' -f1-2 | sort -u)

if [ -z "$changed_plugins" ]; then
  exit 0
fi

marketplace_updated=false

for plugin_dir in $changed_plugins; do
  plugin_json="$plugin_dir/.claude-plugin/plugin.json"

  if [ ! -f "$plugin_json" ]; then
    continue
  fi

  # Get current version
  current_version=$(jq -r '.version' "$plugin_json")

  # Determine bump type from commit message or default to patch
  # Major: breaking changes, removed components
  # Minor: new skills, agents, commands added
  # Patch: bug fixes, documentation updates

  # For now, default to patch bump
  new_version=$(echo "$current_version" | awk -F. '{$NF++; print}' OFS=.)

  # Update plugin.json
  jq --arg version "$new_version" '.version = $version' "$plugin_json" > "$plugin_json.tmp"
  mv "$plugin_json.tmp" "$plugin_json"

  echo "Bumped $plugin_dir version: $current_version → $new_version"

  # Stage the updated plugin.json
  git add "$plugin_json"

  marketplace_updated=true
done

# Bump marketplace version if any plugins changed
if [ "$marketplace_updated" = true ]; then
  marketplace_json=".claude-plugin/marketplace.json"
  current_marketplace=$(jq -r '.metadata.version' "$marketplace_json")
  new_marketplace=$(echo "$current_marketplace" | awk -F. '{$NF++; print}' OFS=.)

  jq --arg version "$new_marketplace" '.metadata.version = $version' "$marketplace_json" > "$marketplace_json.tmp"
  mv "$marketplace_json.tmp" "$marketplace_json"

  echo "Bumped marketplace version: $current_marketplace → $new_marketplace"
  git add "$marketplace_json"
fi

echo "✅ Versions updated automatically"
```

**Configuration in `.claude-plugin/marketplace.json`:**

```json
{
  "metadata": {
    "version": "1.0.0",
    "versioningPolicy": "automatic-patch-bump"
  }
}
```

### Implementation Priority

| Priority | Item                       | Effort | Impact | Dependencies                   |
| -------- | -------------------------- | ------ | ------ | ------------------------------ |
| **P0**   | Pre-commit validation hook | Small  | High   | validation scripts exist       |
| **P0**   | Automated version bumping  | Small  | High   | pre-commit hook                |
| **P1**   | `/validate-plugin` skill   | Medium | High   | plugin-assessor agent          |
| **P1**   | `/create-plugin` skill     | Medium | High   | existing scripts, RT-ICA skill |
| **P2**   | `/extend-plugin` skill     | Medium | Medium | skill-creator, agent-creator   |

### Integration with Methodology

All workflows should follow the process documented in the [SAM methodology](https://github.com/bitflight-devops/stateless-agent-methodology):

1. **RT-ICA Phase** - Reverse Thinking - Information Completeness Assessment

   - Load `rt-ica` skill before planning
   - Identify missing requirements before proceeding
   - Generate comprehensive context manifest

2. **Stateless Agent Methodology (SAM)**

   - Each phase creates artifacts (files) for next phase
   - No reliance on conversation memory
   - Clear handoff protocols between agents

3. **Quality Gates**
   - Each phase has validation checkpoints
   - Cannot proceed to next phase without passing gates
   - Automated validation integrated into workflow

### Validation Improvements Needed

**Current validation tools consolidation:**

Currently validation is spread across:

- `plugin_validator.py` (Python — frontmatter, complexity, links)
- `validate-task-file.sh` (Bash)

**Recommendation:** Create unified `lint-claude-plugin.py`:

- Single Python 3.11+ script
- Cross-platform (Windows/Linux/macOS)
- Pre-commit compatible
- Token-based complexity metrics (not lines)
- Comprehensive validation:
  - Frontmatter schema
  - Plugin.json structure
  - Skill complexity (tokens)
  - Internal link validity
  - Progressive disclosure
  - Tool format correctness
  - Cross-reference validation

---

## Related Plugins

**Optional External Skill Dependency:**

- `prompt-optimization-claude-45` - Used by contextual-ai-documentation-optimizer agent (included in this plugin as agent, requires skill from separate plugin)

**Complementary Plugins:**

- `holistic-linting` - Code quality and linting
- `prompt-optimization-claude-45` - AI-facing documentation optimization
- `python3-development` - Python-specific development patterns

---

## Version History

- **2.6.0** - Consolidated external agents: subagent-refactorer, contextual-ai-documentation-optimizer, plugin-assessor now included
- **2.5.0** - Added skill-creator and write-frontmatter-description skills
- **2.3.0** - Added claude-plugins-reference-2026 and claude-hooks-reference-2026 reference skills
- **2.2.0** - Added claude-skills-overview-2026 reference skill
- **2.1.0** - Added agent-creator skill for creating agents
- **2.0.0** - Merged plugin-creator and plugin-refactor
- **1.0.0** - Initial plugin-creator release

---

## Sources

- Plugin.json: `plugins/plugin-creator/.claude-plugin/plugin.json`
- Skills verified: All 14 SKILL.md files read
- Agents verified: All 6 agent .md files read
- Commands verified: none (commands array empty)
- Scripts verified: All 5 scripts examined
- Official docs: claude-plugins-reference-2026 skill
