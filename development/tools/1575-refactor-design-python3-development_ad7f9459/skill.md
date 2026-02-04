# Refactoring Design Map: python3-development

## Overview

The python3-development plugin requires significant refactoring to address a 1,318-line monolithic SKILL.md covering 5+ distinct domains, non-standard command frontmatter, and orphaned planning documentation. This design splits the monolithic skill into 5 focused skills under 400 lines each, standardizes command schemas, and resolves structural issues while preserving all existing functionality.

## Source Assessment

- **Plugin**: `./plugins/python3-development`
- **Overall Score**: 62/100 (from Phase 1 assessment)
- **Total Refactoring Targets**: 12
  - 1 skill split (into 5 new skills)
  - 6 command frontmatter fixes
  - 1 orphan resolution
  - 4 documentation improvements

---

## Skill Splits

### python3-development Split Plan

**Source**: `./plugins/python3-development/skills/python3-development/SKILL.md`
**Lines**: 1,318
**Domains Identified**:

1. **Orchestration and Workflow** (~200 lines) - Agent delegation, pre-delegation protocol, workflow patterns
2. **Type Safety** (~300 lines) - Generics, Protocols, TypedDict, type narrowing, mypy configuration
3. **CLI Development** (~200 lines) - Typer/Rich patterns, Rich width handling, emoji usage, exception handling
4. **Quality and Linting** (~250 lines) - Linting discovery, format-first workflow, type checker discovery, quality gates
5. **Project Standards** (~200 lines) - Project structure, script trade-offs, PEP 723, asset templates, command usage

**Proposed Split**:

| New Skill          | Scope                                            | Estimated Lines | Sections from Original                       |
| ------------------ | ------------------------------------------------ | --------------- | -------------------------------------------- |
| `python3-workflow` | Agent orchestration, TDD, delegation patterns    | ~200            | Lines 713-846, 924-946, 1227-1246, 1305-1318 |
| `python3-typing`   | Type hints, Generics, Protocols, TypedDict, mypy | ~350            | Lines 274-712                                |
| `python3-cli`      | Typer, Rich, CLI patterns, exception handling    | ~200            | Lines 127-273                                |
| `python3-quality`  | Linting, formatting, quality gates, pre-commit   | ~250            | Lines 941-1129                               |
| `python3-project`  | Project structure, PEP 723, templates, standards | ~200            | Lines 98-126, 1130-1226, 1248-1304           |

### Skill 1: python3-workflow

**Purpose**: Agent orchestration and development workflow patterns for Python projects

**Frontmatter**:

```yaml
---
name: python3-workflow
description: 'Use when: 1. Delegating Python development tasks to agents. 2. Planning TDD, feature addition, or refactoring workflows. 3. Coordinating multi-agent Python implementations. 4. Reviewing orchestration guide requirements. Provides: Agent selection criteria, pre-delegation checklists, workflow patterns (TDD/feature/refactor/debug/review), quality gate sequencing.'
---
```

**Content Structure**:

```markdown
# Python Development Workflow Orchestration

## Role Identification (Mandatory)
[Keep existing role identification section - lines 11-35]

<section ROLE_TYPE="orchestrator">

## Agent Orchestration (Orchestrator Only)

### STOP - MANDATORY PRE-DELEGATION PROTOCOL
[Lines 719-731 - Pre-delegation checklist]

### Why This Protocol Exists
[Lines 733-740]

### DO NOT Pre-Gather Data for Agents
[Lines 754-776]

### Delegation Pattern
[Lines 778-805 - Delegation table]

### Required Reading (MANDATORY)
[Lines 806-819]

### Quick Reference Example
[Lines 820-846]

## Core Workflows (Orchestrator Only)
[Lines 924-940]

## Common Patterns to Follow (Orchestrator Only)
[Lines 1227-1246]

</section>

## Summary for All Roles
[Brief summary with link to orchestration guide reference]
```

**Shared References**:

- `./references/python-development-orchestration.md` (primary reference)
- `./references/accessing_online_resources.md` (symlink)

**Migration Notes**:

- This skill focuses ONLY on orchestration patterns
- Sub-agents do NOT need this skill - they receive task instructions directly
- Must include ROLE_TYPE sections to maintain orchestrator-only guidance

---

### Skill 2: python3-typing

**Purpose**: Comprehensive type annotation patterns for Python 3.11+

**Frontmatter**:

```yaml
---
name: python3-typing
description: 'Use when: 1. Adding type hints to Python code. 2. Working with Generics, Protocols, or TypedDict. 3. Configuring mypy or fixing type errors. 4. Using type narrowing or TypeGuard patterns. 5. Choosing between attrs, dataclasses, and pydantic. Provides: Python 3.11+ native type hints, Generic patterns, Protocol structural typing, TypedDict schemas, type narrowing, mypy configuration.'
---
```

**Content Structure**:

```markdown
# Python Type Safety with Mypy

## Type Safety Requirements
[Brief intro - requirement for comprehensive typing]

## When to Use Generics
[Lines 281-351 - Generic patterns with TypeVar]

## When to Use Protocols
[Lines 365-436 - Protocol patterns for structural typing]

## TypedDict for Dictionary Typing
[Lines 438-506 - TypedDict patterns]

## Type Narrowing
[Lines 508-569 - isinstance, None checks, TypeGuard, TypeIs]

## attrs vs dataclasses vs pydantic
[Lines 571-639 - Decision matrix and patterns]

## Additional Mypy Features
[Lines 641-677 - Generic dataclasses, Self type]

## Mypy Configuration Best Practices
[Lines 679-711 - pyproject.toml configuration]
```

**Shared References**:

- `./references/mypy-docs/generics.rst`
- `./references/mypy-docs/protocols.rst`
- `./references/mypy-docs/typed_dict.rst`
- `./references/mypy-docs/type_narrowing.rst`
- `./references/mypy-docs/additional_features.rst`
- `./references/tool-library-registry.md` (for mypy config section)

**Migration Notes**:

- Largest skill due to comprehensive type patterns
- All mypy-docs references remain unchanged
- No ROLE_TYPE sections needed - applies to all roles

---

### Skill 3: python3-cli

**Purpose**: CLI development patterns with Typer and Rich

**Frontmatter**:

```yaml
---
name: python3-cli
description: 'Use when: 1. Building CLI applications with Typer or Rich. 2. Fixing Rich Panel/Table width wrapping issues. 3. Handling exceptions in CLI apps. 4. Using Rich emoji tokens. 5. Creating CLI output for CI/non-TTY environments. Provides: Typer+Rich integration patterns, width handling solutions, exception handling patterns, Rich emoji usage, CI-compatible output.'
---
```

**Content Structure**:

```markdown
# Python CLI Development with Typer and Rich

## Script Dependency Trade-offs
[Lines 127-156 - When to use Typer+Rich vs stdlib]

## Rich Panel and Table Width Handling
[Lines 157-217 - get_rendered_width() pattern]

## Rich Emoji Usage
[Lines 219-247 - Token patterns]

## Python Exception Handling
[Lines 249-271 - Fail-fast principle]

## Exception Handling in CLI Applications
[Brief section pointing to reference]
```

**Shared References**:

- `./references/exception-handling.md` (critical for Typer exception chain prevention)
- `./assets/typer_examples/index.md`
- `./assets/typer_examples/console_no_wrap_example.py`
- `./assets/typer_examples/console_containers_no_wrap.py`

**Migration Notes**:

- Keep executable examples in assets/
- Exception handling reference is critical for Typer apps
- No ROLE_TYPE sections needed - applies to all roles

---

### Skill 4: python3-quality

**Purpose**: Linting, formatting, and quality gate execution

**Frontmatter**:

```yaml
---
name: python3-quality
description: 'Use when: 1. Running linting or formatting on Python code. 2. Discovering project linters (pre-commit, prek, ruff, mypy). 3. Resolving linting errors at root cause. 4. Setting up quality gates. 5. Working with basedpyright, pyright, or mypy type checkers. Provides: Linting discovery protocol, format-first workflow, type checker detection, quality gate sequences, linting exception conditions.'
---
```

**Content Structure**:

```markdown
# Python Code Quality and Linting

## Linting Discovery Protocol
[Lines 944-1046 - Discovery sequence, format-first, type checker discovery]

## Quality Gates
[Lines 1047-1079 - Gate sequence with type checker detection]

## Linting Exception Conditions
[Lines 1087-1129 - When to ignore vs fix]
```

**Shared References**:

- `./references/tool-library-registry.md` (for tool configurations)
- `./references/user-project-conventions.md` (for ruff/mypy settings)

**Migration Notes**:

- Critical integration with holistic-linting skill
- Must detect both pre-commit and prek tools
- Include CI compatibility verification

---

### Skill 5: python3-project

**Purpose**: Project structure, standards, and configuration

**Frontmatter**:

```yaml
---
name: python3-project
description: 'Use when: 1. Creating new Python projects or packages. 2. Setting up pyproject.toml configuration. 3. Using PEP 723 inline script metadata. 4. Copying asset templates (version.py, pre-commit config). 5. Understanding standard project structure. Provides: Project layout standards, pyproject.toml templates, PEP 723 guidance, asset template usage, command usage (/modernpython, /shebangpython).'
---
```

**Content Structure**:

```markdown
# Python Project Standards

## Python Development Standards
[Lines 98-126 - Core concepts, docstring standard, template variables]

## Standard Project Structure
[Lines 1130-1162 - Directory layout, hatchling config]

## Using Asset Templates
[Lines 1182-1226 - Template usage instructions]

## Command Usage
[Lines 848-923 - /modernpython and /shebangpython]

## Integration
[Lines 1163-1181 - External reference example]

## Detailed Documentation
[Lines 1248-1292 - Reference links]
```

**Shared References**:

- `./references/PEP723.md`
- `./references/user-project-conventions.md`
- `./references/tool-library-registry.md` (template variable reference)
- `./references/modern-modules.md`
- `./assets/` directory (all templates)

**Migration Notes**:

- Keep all asset templates in existing location
- Commands reference external files in `~/.claude/commands/`
- Include reference navigation helpers (grep commands)

---

## Shared Reference Strategy

All 5 new skills share access to the existing `references/` directory:

| Reference File                        | Used By Skills                                   |
| ------------------------------------- | ------------------------------------------------ |
| `python-development-orchestration.md` | python3-workflow                                 |
| `exception-handling.md`               | python3-cli                                      |
| `mypy-docs/*.rst` (5 files)           | python3-typing                                   |
| `tool-library-registry.md`            | python3-typing, python3-quality, python3-project |
| `user-project-conventions.md`         | python3-quality, python3-project                 |
| `PEP723.md`                           | python3-project                                  |
| `modern-modules.md`                   | python3-project                                  |
| `modern-modules/*.md` (18 files)      | python3-project (on-demand)                      |
| `api_reference.md`                    | python3-project                                  |
| `accessing_online_resources.md`       | python3-workflow                                 |

**Directory Structure After Split**:

```
plugins/python3-development/
├── skills/
│   ├── python3-workflow/
│   │   └── SKILL.md
│   ├── python3-typing/
│   │   └── SKILL.md
│   ├── python3-cli/
│   │   └── SKILL.md
│   ├── python3-quality/
│   │   └── SKILL.md
│   └── python3-project/
│       ├── SKILL.md
│       ├── references/           # Moved from python3-development
│       │   ├── python-development-orchestration.md
│       │   ├── exception-handling.md
│       │   ├── PEP723.md
│       │   ├── user-project-conventions.md
│       │   ├── tool-library-registry.md
│       │   ├── modern-modules.md
│       │   ├── modern-modules/
│       │   ├── mypy-docs/
│       │   ├── api_reference.md
│       │   └── accessing_online_resources.md -> symlink
│       └── assets/               # Moved from python3-development
│           ├── typer_examples/
│           ├── version.py
│           ├── hatch_build.py
│           └── ...
├── commands/
│   ├── development/
│   └── testing/
├── plugin.json
└── README.md
```

**Rationale for References Location**:

- References stay with `python3-project` as the primary skill
- Other skills reference files using relative paths: `../python3-project/references/`
- Alternative: Create shared `references/` at plugin level and update all paths

---

## Command Frontmatter Fixes

### Current Schema Issues

All 6 commands use non-standard frontmatter:

| Field          | Current (Non-Standard) | Required (Standard)     |
| -------------- | ---------------------- | ----------------------- |
| `title`        | Present                | Remove (not in schema)  |
| `command_type` | Present                | Remove (not in schema)  |
| `description`  | Missing                | Required for /help menu |

### Commands to Fix

#### 1. create-feature-task.md

**Source**: `./plugins/python3-development/commands/development/create-feature-task.md`

**Current Frontmatter**:

```yaml
---
title: "Create Feature Development Task"
description: "Set up comprehensive feature development task with proper tracking"
command_type: "development"
last_updated: "2025-11-02"
related_docs:
  - "./use-command-template.md"
  - "../../references/python-development-orchestration.md"
---
```

**Target Frontmatter**:

```yaml
---
name: create-feature-task
description: 'Create structured feature development task with tracking. Use when: starting new feature work, setting up task phases, creating development checkpoints.'
argument-hint: "[feature-name]"
---
```

**Changes**:

- Remove `title` (use `name` or let it default from filename)
- Remove `command_type` (not in schema)
- Remove `last_updated` (not in schema)
- Remove `related_docs` (put in content as markdown links)
- Add `argument-hint` for better autocomplete
- Rewrite description with trigger keywords

---

#### 2. use-command-template.md

**Source**: `./plugins/python3-development/commands/development/use-command-template.md`

**Target Frontmatter**:

```yaml
---
name: use-command-template
description: 'Create new slash command from template. Use when: adding custom commands, creating development shortcuts, setting up workflow automation.'
argument-hint: "[command-name]"
---
```

---

#### 3. command-template.md (template file)

**Source**: `./plugins/python3-development/commands/development/templates/command-template.md`

**Target Frontmatter**:

```yaml
---
name: command-template
description: 'Template for creating new slash commands. Not for direct invocation - copy and customize for new commands.'
user-invocable: false
---
```

**Changes**:

- Add `user-invocable: false` since it's a template, not a command

---

#### 4. analyze-test-failures.md

**Source**: `./plugins/python3-development/commands/testing/analyze-test-failures.md`

**Target Frontmatter**:

```yaml
---
name: analyze-test-failures
description: 'Analyze failing tests with investigative approach. Use when: tests fail unexpectedly, determining if test or implementation is wrong, debugging test assertions.'
argument-hint: "[test-file-or-function]"
---
```

---

#### 5. comprehensive-test-review.md

**Source**: `./plugins/python3-development/commands/testing/comprehensive-test-review.md`

**Target Frontmatter**:

```yaml
---
name: comprehensive-test-review
description: 'Review test suite for completeness and quality. Use when: auditing test coverage, reviewing test patterns, assessing test suite health.'
argument-hint: "[test-directory]"
---
```

---

#### 6. test-failure-mindset.md

**Source**: `./plugins/python3-development/commands/testing/test-failure-mindset.md`

**Target Frontmatter**:

```yaml
---
name: test-failure-mindset
description: 'Load test debugging mindset guidelines. Use when: approaching test failures, avoiding assumption that tests are wrong, systematic failure analysis.'
---
```

---

## Orphan Resolution

### planning/reference-document-architecture.md

**Source**: `./plugins/python3-development/skills/python3-development/planning/reference-document-architecture.md`

**Classification**: Historical Proposal (Working Notes)

**Current Status**:

- Not linked from SKILL.md or any reference file
- Frontmatter explicitly states: `status: "historical-proposal"`, `implementation_status: "not-implemented"`
- Contains architectural proposal that was never implemented

**Resolution**: Archive or Remove

**Recommendation**: Move to documentation or remove entirely

**Option A - Archive** (Recommended):

1. Move to `./plugins/python3-development/docs/archive/reference-document-architecture.md`
2. Add link from README.md under "Historical Documentation" section

**Option B - Remove**:

1. Delete the file since it's marked as not-implemented
2. The actual skill structure differs significantly from this proposal

**Implementation**:

```bash
# Option A - Archive
mkdir -p ./plugins/python3-development/docs/archive/
mv ./plugins/python3-development/skills/python3-development/planning/reference-document-architecture.md \
   ./plugins/python3-development/docs/archive/
rmdir ./plugins/python3-development/skills/python3-development/planning/

# Update README.md to mention archived docs
```

---

## Documentation Improvements

### 1. plugin.json Description

**Target**: `./plugins/python3-development/plugin.json`

**Current Issue**: File does not exist (404 error in assessment read)

**Action**: Create plugin.json with proper metadata

**Target Content**:

```json
{
  "name": "python3-development",
  "version": "2.0.0",
  "description": "Python 3.11+ development patterns, type safety, CLI tools, and quality workflows",
  "author": "User",
  "license": "MIT",
  "keywords": ["python", "typing", "typer", "rich", "mypy", "ruff", "pytest"],
  "skills": [
    "python3-workflow",
    "python3-typing",
    "python3-cli",
    "python3-quality",
    "python3-project"
  ],
  "commands": {
    "development": ["create-feature-task", "use-command-template"],
    "testing": ["analyze-test-failures", "comprehensive-test-review", "test-failure-mindset"]
  }
}
```

---

### 2. Reference Files Back-Links

**Target**: All 26 reference files in `./references/`

**Current Issue**: No navigation headers pointing back to parent skill

**Action**: Add navigation header to each reference file

**Template**:

```markdown
---
[existing frontmatter]
parent_skill: python3-project
---

> **Navigation**: [Back to Python Project Standards](../SKILL.md) | [Reference Index](./modern-modules.md)

[existing content]
```

**Priority**: Low - Can be done incrementally

---

### 3. Skill Description Optimization

All 5 new skills need trigger-optimized descriptions (already specified in Skill Splits section above).

**Verification Checklist**:

- [ ] Under 1024 characters
- [ ] Contains "Use when:" trigger phrases
- [ ] Contains "Provides:" capability list
- [ ] No YAML multiline indicators
- [ ] Single-quoted if contains special characters

---

## Dependency Map

```
                    ┌─────────────────────────────────────────────────────┐
                    │                  PHASE 1: STRUCTURE                  │
                    │   (Create directories, plugin.json, archive orphan)  │
                    └─────────────────────────────────────────────────────┘
                                             │
                                             ▼
         ┌───────────────────────────────────┴───────────────────────────────────┐
         │                                                                       │
         ▼                                                                       ▼
┌─────────────────────────┐                                       ┌─────────────────────────┐
│   PHASE 2A: SKILL SPLIT │                                       │  PHASE 2B: CMD FIXES    │
│   (5 new SKILL.md files)│                                       │  (6 frontmatter fixes)  │
└─────────────────────────┘                                       └─────────────────────────┘
         │                                                                       │
         │  ┌─────────────┬─────────────┬─────────────┬─────────────┐           │
         │  │             │             │             │             │           │
         ▼  ▼             ▼             ▼             ▼             ▼           │
    python3-      python3-      python3-      python3-      python3-            │
    workflow      typing        cli           quality       project             │
         │             │             │             │             │              │
         └─────────────┴─────────────┴─────────────┴─────────────┴──────────────┘
                                             │
                                             ▼
                    ┌─────────────────────────────────────────────────────┐
                    │                  PHASE 3: REFERENCES                 │
                    │      (Move references to python3-project, update     │
                    │       relative paths in other skills)                │
                    └─────────────────────────────────────────────────────┘
                                             │
                                             ▼
                    ┌─────────────────────────────────────────────────────┐
                    │                  PHASE 4: VALIDATION                 │
                    │      (Test skill loading, verify cross-references,   │
                    │       run pre-commit, update README)                 │
                    └─────────────────────────────────────────────────────┘
```

---

## Parallelization Opportunities

### Phase 2A: Skill Split (Parallelizable)

All 5 skills can be created simultaneously as they have no dependencies on each other during creation:

**Group A - No dependencies between tasks**:

| Task                             | Agent            | Target File                        |
| -------------------------------- | ---------------- | ---------------------------------- |
| Create python3-workflow SKILL.md | skill-refactorer | `skills/python3-workflow/SKILL.md` |
| Create python3-typing SKILL.md   | skill-refactorer | `skills/python3-typing/SKILL.md`   |
| Create python3-cli SKILL.md      | skill-refactorer | `skills/python3-cli/SKILL.md`      |
| Create python3-quality SKILL.md  | skill-refactorer | `skills/python3-quality/SKILL.md`  |
| Create python3-project SKILL.md  | skill-refactorer | `skills/python3-project/SKILL.md`  |

### Phase 2B: Command Fixes (Parallelizable)

All 6 commands can be fixed simultaneously:

**Group B - No dependencies between tasks**:

| Task                                      | Agent                    | Target File                                          |
| ----------------------------------------- | ------------------------ | ---------------------------------------------------- |
| Fix create-feature-task frontmatter       | claude-context-optimizer | `commands/development/create-feature-task.md`        |
| Fix use-command-template frontmatter      | claude-context-optimizer | `commands/development/use-command-template.md`       |
| Fix command-template frontmatter          | claude-context-optimizer | `commands/development/templates/command-template.md` |
| Fix analyze-test-failures frontmatter     | claude-context-optimizer | `commands/testing/analyze-test-failures.md`          |
| Fix comprehensive-test-review frontmatter | claude-context-optimizer | `commands/testing/comprehensive-test-review.md`      |
| Fix test-failure-mindset frontmatter      | claude-context-optimizer | `commands/testing/test-failure-mindset.md`           |

### Sequential Dependencies

- Phase 1 MUST complete before Phase 2 (directory structure required)
- Phase 2 MUST complete before Phase 3 (skills must exist before moving references)
- Phase 3 MUST complete before Phase 4 (paths must be correct before validation)

---

## Implementation Task List

### Phase 1: Structure Setup

| ID  | Task                                                  | Priority | Agent        | Estimated Time |
| --- | ----------------------------------------------------- | -------- | ------------ | -------------- |
| 1.1 | Create 5 skill directories under `skills/`            | Critical | orchestrator | 5 min          |
| 1.2 | Create `plugin.json` with metadata                    | Critical | orchestrator | 10 min         |
| 1.3 | Archive `planning/reference-document-architecture.md` | Medium   | orchestrator | 5 min          |
| 1.4 | Create `docs/archive/` directory                      | Medium   | orchestrator | 2 min          |

### Phase 2A: Skill Split

| ID  | Task                               | Priority | Agent            | Estimated Time |
| --- | ---------------------------------- | -------- | ---------------- | -------------- |
| 2.1 | Create `python3-workflow/SKILL.md` | Critical | skill-refactorer | 30 min         |
| 2.2 | Create `python3-typing/SKILL.md`   | Critical | skill-refactorer | 45 min         |
| 2.3 | Create `python3-cli/SKILL.md`      | Critical | skill-refactorer | 25 min         |
| 2.4 | Create `python3-quality/SKILL.md`  | Critical | skill-refactorer | 30 min         |
| 2.5 | Create `python3-project/SKILL.md`  | Critical | skill-refactorer | 35 min         |

### Phase 2B: Command Fixes

| ID  | Task                            | Priority | Agent                    | Estimated Time |
| --- | ------------------------------- | -------- | ------------------------ | -------------- |
| 2.6 | Fix 6 command frontmatter files | High     | claude-context-optimizer | 20 min         |

### Phase 3: Reference Migration

| ID  | Task                                              | Priority | Agent            | Estimated Time |
| --- | ------------------------------------------------- | -------- | ---------------- | -------------- |
| 3.1 | Move `references/` to `python3-project/`          | Critical | orchestrator     | 10 min         |
| 3.2 | Move `assets/` to `python3-project/`              | Critical | orchestrator     | 5 min          |
| 3.3 | Update relative paths in 4 other skills           | Critical | skill-refactorer | 20 min         |
| 3.4 | Delete old `python3-development/` skill directory | Critical | orchestrator     | 2 min          |

### Phase 4: Validation

| ID  | Task                                        | Priority | Agent              | Estimated Time |
| --- | ------------------------------------------- | -------- | ------------------ | -------------- |
| 4.1 | Test plugin with `--plugin-dir`             | Critical | orchestrator       | 2 min          |
| 4.2 | Test skill loading with `/skill-name`       | Critical | orchestrator       | 15 min         |
| 4.3 | Verify all reference file links resolve     | High     | orchestrator       | 10 min         |
| 4.4 | Run pre-commit on all modified files        | High     | orchestrator       | 5 min          |
| 4.5 | Update `README.md` with new skill structure | Medium   | plugin-docs-writer | 20 min         |

---

## Risk Assessment

| Risk                                       | Likelihood | Impact | Mitigation                                                   |
| ------------------------------------------ | ---------- | ------ | ------------------------------------------------------------ |
| Broken reference links after migration     | High       | Medium | Validate all links in Phase 4 before deleting old skill      |
| Skill truncation from `<available_skills>` | Medium     | Medium | Keep descriptions under 300 chars, test with 5 skills active |
| Missing content during split               | Medium     | High   | Use checklist to verify all lines accounted for              |
| Plugin loading issues                      | Low        | High   | Test skill loading immediately with `--plugin-dir`           |

---

## Success Criteria

### Quantitative

- [ ] Each new skill is under 400 lines (target: 200-350)
- [ ] All 5 new skills load without errors
- [ ] All 6 commands appear in `/help` menu with descriptions
- [ ] All reference file links resolve correctly
- [ ] Pre-commit passes on all modified files
- [ ] No content lost from original SKILL.md

### Qualitative

- [ ] Each skill has single, focused domain
- [ ] Skill descriptions contain clear trigger keywords
- [ ] Progressive disclosure pattern maintained (SKILL.md -> references/)
- [ ] Orchestrator-only sections preserved in python3-workflow
- [ ] Asset templates remain accessible from python3-project

---

## Rollback Plan

If issues discovered after migration:

1. **Restore from git**: All changes tracked in git, can revert to pre-refactor state
2. **Keep old skill temporarily**: Rename to `python3-development-legacy/` during transition
3. **Symlink fallback**: Update `~/.claude/skills/` symlink to point to legacy skill

---

## Document Metadata

- **Created**: 2026-01-22
- **Author**: python-cli-design-spec agent
- **Assessment Source**: Plugin assessment report (Phase 1)
- **Target Plugin**: `./plugins/python3-development`
- **Target Score**: 85+/100 (from 62/100)
