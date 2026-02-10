# Refactoring Design Map: holistic-linting

## Overview

The holistic-linting plugin has a single 1221-line SKILL.md covering 4 distinct domains (orchestrator workflows, sub-agent procedures, linter-specific resolution workflows, and bundled resources). This refactoring splits SKILL.md into focused, role-appropriate skills while optimizing both agents and fixing documentation issues. The plugin's reference file organization is exemplary and requires no changes.

## Source Assessment

- Plugin: ./plugins/holistic-linting
- Overall Score: 68/100 (penalized for oversized SKILL.md and agent model choice)
- Total Refactoring Targets: 6

## Corrected Assessment Notes

The Phase 1 assessment incorrectly flagged `./scripts/` as a broken link. The scripts directory exists at `skills/holistic-linting/scripts/` and the `./scripts/` link in SKILL.md resolves correctly relative to SKILL.md's location. However, lines 1005-1011 use `python holistic-linting/scripts/install-agents.py` which is relative to the skills directory parent, not to SKILL.md — this is a documentation accuracy issue (wrong relative path in code examples), not a broken file link.

## Skill Split

### holistic-linting SKILL.md Split Plan

**Source**: ./plugins/holistic-linting/skills/holistic-linting/SKILL.md
**Lines**: 1221
**Domains Identified**: 4

1. **Core skill (router/overview)** — Lines 1-48, 247-375, 993-1222
2. **Orchestrator delegation guide** — Lines 49-245 (inside `<section ROLE_TYPE="orchestrator">`)
3. **Linter resolution workflows** — Lines 376-991 (Ruff, Mypy, Pyright procedures)

**Proposed Split**:

| New Skill | Scope | Source Lines | Est. Lines |
|-----------|-------|-------------|------------|
| holistic-linting (trimmed) | Router: purpose, when to use, linter detection, running formatters, bundled resources, examples, best practices | 1-48, 247-375, 993-1222 | ~430 |
| holistic-linting-orchestrator | Orchestrator-only delegation workflow, anti-patterns, report reading | 49-245 | ~200 |
| holistic-linting-resolver | Sub-agent resolution workflows for Ruff, Mypy, Pyright including suppression gates | 376-991 | ~616 → target ~500 via deduplication |

**Deduplication Opportunity in Resolver Skill**:

The three resolution workflows (Ruff lines 380-513, Mypy lines 514-707, Pyright lines 709-957) share identical patterns:
- Suppression Gate section (duplicated 3x, ~12 lines each = 36 lines)
- "Load python3-development skill" step (duplicated 3x)
- "Check Architectural Context" step (similar across all 3)
- Verification step pattern (similar across all 3)

Extract shared methodology as a "Common Resolution Steps" section (~40 lines), then each linter-specific section only documents what's unique. This should reduce ~616 lines to ~500 lines.

**Shared References**: All three skills share the same `references/` directory. The existing progressive disclosure hierarchy (SKILL.md → index files → rule files) remains unchanged. Only SKILL.md links to the index files; the two new skills do not need direct references to the rules knowledge base since:
- `holistic-linting-orchestrator` delegates to agents (doesn't need rules)
- `holistic-linting-resolver` will reference the rules knowledge base through the main skill's established paths

**Migration Notes**:
- The main `holistic-linting` skill is referenced by both agents (linting-root-cause-resolver loads it via `Skill(command: "holistic-linting")`). After split, this agent should load `holistic-linting` (which loads resolver content or references it).
- The `<section ROLE_TYPE="orchestrator">` XML tag on line 49 already delineates the orchestrator content cleanly.
- The `/lint` command (commands/lint.md) references the main skill workflow — no changes needed since the main skill remains the entry point.

**Skill Directory Structure After Split**:

```
plugins/holistic-linting/
├── .claude-plugin/plugin.json          # Updated to list 3 skills
├── commands/lint.md                     # Unchanged
├── agents/
│   ├── linting-root-cause-resolver.md   # Updated description
│   └── post-linting-architecture-reviewer.md  # Model + description updated
├── skills/
│   ├── holistic-linting/
│   │   ├── SKILL.md                     # Trimmed to ~430 lines (router)
│   │   ├── scripts/                     # Unchanged
│   │   └── references/                  # Unchanged (all rules, mypy-docs)
│   ├── holistic-linting-orchestrator/
│   │   └── SKILL.md                     # ~200 lines (orchestrator delegation)
│   └── holistic-linting-resolver/
│       └── SKILL.md                     # ~500 lines (resolution workflows)
└── README.md                            # Updated
```

## Agent Optimizations

### post-linting-architecture-reviewer Optimization

**Source**: ./plugins/holistic-linting/agents/post-linting-architecture-reviewer.md
**Lines**: 182

**Current Description**: "Architectural review after linting-root-cause-resolver completes. Verifies resolution quality, examines artifacts in .claude/reports/, checks fixes align with codebase patterns, and identifies systemic improvements. Trigger after linting resolution."

**Issues**:
1. Model set to `haiku` — inappropriate for architectural reasoning (reasoning-heavy task)
2. Description lacks keywords for design principles and type safety analysis
3. Lines 116-165 contain a markdown code fence that is never properly closed (the closing ``` on line 165 is inside a nested code block context, and line 182 has a stray ```)

**Proposed Changes**:

1. **Change model from `haiku` to `inherit`** (line 4)
   - Architectural review requires reliable reasoning about design patterns
   - `inherit` uses the orchestrator's model, appropriate for reasoning tasks

2. **Enhanced description**: "Architectural review after linting-root-cause-resolver completes. Verifies resolution quality, examines artifacts in .claude/reports/, checks fixes align with codebase patterns and design principles, validates type safety improvements, and identifies systemic improvements. Use after linting resolution to assess SOLID compliance, code organization, and broader architectural impact."

3. **Fix markdown structure**: The agent file has a structural issue — the code fence block starting at line 116 (````markdown) isn't cleanly closed. Lines 164-165 and 182 have stray fence markers. Clean up the closing structure.

### linting-root-cause-resolver — No Changes Needed

**Source**: ./plugins/holistic-linting/agents/linting-root-cause-resolver.md
**Lines**: 163
**Assessment**: 10/10 description quality, correct model (`inherit`), clear instructions
**Action**: None required

## Documentation Improvements

### 1. Fix install-agents.py Path in SKILL.md

**Source**: ./plugins/holistic-linting/skills/holistic-linting/SKILL.md lines 1003-1012
**Current**: `python holistic-linting/scripts/install-agents.py --scope user`
**Issue**: Path is relative to skills directory parent, inconsistent with other script references that use `./scripts/`
**Fix**: Change to `python ./scripts/install-agents.py --scope user` (or `uv run ./scripts/install-agents.py --scope user` for consistency)

### 2. plugin.json Metadata Enhancement

**Source**: ./plugins/holistic-linting/.claude-plugin/plugin.json
**Current**: Missing keywords, repository, license, homepage
**Proposed additions**:

```json
{
  "keywords": ["linting", "code-quality", "ruff", "mypy", "bandit", "python", "formatting", "type-checking"],
  "repository": "https://github.com/bitflight-devops/claude_skills/tree/main/plugins/holistic-linting",
  "license": "MIT"
}
```

### 3. Update plugin.json skills array

After the skill split, update `skills` in plugin.json:

```json
{
  "skills": [
    "./skills/holistic-linting",
    "./skills/holistic-linting-orchestrator",
    "./skills/holistic-linting-resolver"
  ]
}
```

## Orphan Resolution

No orphaned files identified. All 52 reference files are properly indexed through the hierarchical index system. The mypy-docs/*.rst files are referenced in SKILL.md lines 526-529 as cached documentation sources.

## Dependency Map

```
Task 1: Split SKILL.md into 3 skills
  ├── No dependencies (can start immediately)
  └── Blocks: Task 3 (plugin.json update), Task V1, Task V2

Task 2: Optimize post-linting-architecture-reviewer agent
  ├── No dependencies (can start immediately)
  └── Blocks: Task V1

Task 3: Update plugin.json (metadata + skills array)
  ├── Depends on: Task 1 (needs new skill directories to exist)
  └── Blocks: Task V1

Task 4: Fix install-agents.py paths in SKILL.md
  ├── Depends on: Task 1 (will be part of SKILL.md content during split)
  └── Can be merged into Task 1

Task V1: Validate plugin structure
  ├── Depends on: Tasks 1, 2, 3
  └── Blocks: Task V2

Task V2: Update plugin documentation (README.md)
  ├── Depends on: Task V1
  └── No blockers
```

## Parallelization Opportunities

Tasks that can run simultaneously:

- **Group A**: Task 1 (Skill Split) + Task 2 (Agent Optimization) — No shared files
- **Group B**: Task 3 (plugin.json update) — Depends on Task 1 completing
- **Group C**: Task V1 (Validation) — Depends on Tasks 1, 2, 3
- **Group D**: Task V2 (README update) — Depends on Task V1
