# Plugin Refactoring Plan Index

This index tracks all plugin refactoring projects, linking to detailed task breakdowns and tracking progress from initial assessment through completion.

---

## Active Refactoring Projects

| Plugin              | Task File                                                                        | Status         | Score Before | Score After | Phase |
| ------------------- | -------------------------------------------------------------------------------- | -------------- | ------------ | ----------- | ----- |
| python3-development | [tasks-refactor-python3-development.md](./tasks-refactor-python3-development.md) | ❌ NOT STARTED | 62/100       | TBD         | 3     |

---

## Completed Refactoring Projects

| Plugin | Task File | Completion Date | Score Improvement | Final Score |
| ------ | --------- | --------------- | ----------------- | ----------- |
| _None_ | _-_       | _-_             | _-_               | _-_         |

---

## Refactoring Workflow

### Phase 1: Assessment

**Agent**: `plugin-assessor`
**Output**: Assessment report identifying refactoring targets

**Triggers refactoring when**:

- Overall score <70/100
- Monolithic skills >500 lines
- Non-standard frontmatter detected
- Orphaned files found
- Missing critical documentation

### Phase 2: Design

**Agent**: `python-cli-design-spec`
**Output**: Refactoring design document with skill splits, frontmatter fixes, orphan resolution

**Key Deliverables**:

- Skill split plan with line number mappings
- Frontmatter standardization specifications
- Orphan resolution strategy
- Dependency map
- Parallelization opportunities

### Phase 3: Task Breakdown (This Phase)

**Agent**: `task-planner`
**Output**: Detailed task file with dependency-based ordering

**Key Deliverables**:

- Task file with acceptance criteria
- Dependency mapping
- Parallelization identification
- Agent assignment per task
- Verification steps

### Phase 4: Execution

**Agent**: `orchestrator` (delegates to specialist agents)
**Output**: Refactored plugin structure

**Specialist Agents**:

- `skill-refactorer` - Splits monolithic skills
- `subagent-refactorer` - Optimizes agent files
- `claude-context-optimizer` - Improves AI-facing docs
- `plugin-docs-writer` - Generates user documentation
- `plugin-assessor` - Validates final structure

### Phase 5: Validation

**Agent**: `plugin-assessor` + `orchestrator`
**Output**: Updated assessment, passing pre-commit, working symlinks

**Validation Checklist**:

- [ ] Plugin structure passes plugin-assessor
- [ ] All skills load via symlinks
- [ ] Commands appear in `/help` menu
- [ ] Pre-commit hooks pass
- [ ] Score improved from baseline
- [ ] No orphaned files remain

---

## Task File Template

Each plugin refactoring task file follows this structure:

### Required Sections

1. **Task Definitions** - Each with:

   - Status (❌ NOT STARTED | 🔄 IN PROGRESS | ✅ COMPLETE)
   - Dependencies (explicit task IDs or "None")
   - Priority (1-5 integer)
   - Complexity (Low/Medium/High)
   - Agent (which specialist agent executes)
   - Acceptance Criteria (minimum 3 specific criteria)
   - Required Inputs (design spec sections, source files)
   - Expected Outputs (file paths created/modified)
   - Parallelization markers (which tasks can run concurrently)
   - Verification Steps (minimum 3 testable steps)

2. **Dependency Summary** - Phases with parallelization analysis

3. **Success Metrics** - Quantitative and qualitative targets

### Task Naming Convention

- **Phase 1 (Foundation)**: Tasks 1.1, 1.2, 1.3...
- **Phase 2A (Skill Splits)**: Tasks 2.1, 2.2, 2.3...
- **Phase 2B (Command Fixes)**: Task 2.6
- **Phase 3 (Migration)**: Tasks 3.1, 3.2, 3.3, 3.4
- **Phase 4 (Validation)**: Tasks V1, V2, V3, V4

---

## Progress Tracking Protocol

### Starting a Refactoring Project

1. Plugin assessment completed (Phase 1)
2. Refactoring design document created (Phase 2)
3. Task breakdown file created (Phase 3) ← **Current step**
4. Add entry to "Active Refactoring Projects" table above
5. Link to task file with relative path

### Updating Task Status

**Rule**: Update task file status in real-time as work progresses

```markdown
**Status**: ❌ NOT STARTED → 🔄 IN PROGRESS → ✅ COMPLETE
```

**Update Triggers**:

- ❌→🔄: When agent begins work on task
- 🔄→✅: When all acceptance criteria met and verification steps pass
- Never mark ✅ if errors/blockers remain

### Completing a Refactoring Project

1. All tasks marked ✅ COMPLETE
2. Final plugin assessment shows score improvement
3. Move entry from "Active" to "Completed" table
4. Record completion date and final score
5. Calculate score improvement (Final - Before)

---

## Refactoring Anti-Patterns to Avoid

### ❌ Temporal Language in Tasks

```markdown
# BAD
Week 1: Create skill directories
Sprint 1: Split monolithic skill (2 weeks)
```

```markdown
# GOOD
**Dependencies**: None (can run in parallel with Task 1.2, 1.3)
**Priority**: 1 (foundational - blocks Phase 2)
```

### ❌ Vague Acceptance Criteria

```markdown
# BAD
Acceptance Criteria:
1. Skill is well-structured
2. File is readable
```

```markdown
# GOOD
Acceptance Criteria:
1. File length under 400 lines (target: ~350 lines)
2. Frontmatter description contains "Use when:" and "Provides:" sections
3. All mypy-docs references preserved with correct relative paths
```

### ❌ Missing Parallelization Analysis

```markdown
# BAD
**Can Parallelize With**: Other tasks

# GOOD
**Can Parallelize With**: Task 2.2, Task 2.3, Task 2.4, Task 2.5
**Reason**: Each skill writes to different file - no shared dependencies during creation
```

### ❌ Incomplete Verification Steps

```markdown
# BAD
Verification Steps:
1. Check if file exists

# GOOD
Verification Steps:
1. Run `wc -l ./plugins/python3-development/skills/python3-typing/SKILL.md` - should be under 400 lines
2. Run `grep -c "mypy-docs" ./plugins/python3-development/skills/python3-typing/SKILL.md` - should find multiple references
3. Verify all sections from design spec present: Generics, Protocols, TypedDict, Type Narrowing
```

---

## Metrics and Reporting

### Score Improvement Targets

- **Baseline <50/100**: Target +40 points (reach 90+)
- **Baseline 50-69/100**: Target +20-30 points (reach 80+)
- **Baseline 70-79/100**: Target +10-15 points (reach 85+)
- **Baseline 80+/100**: Target +5-10 points (reach 90+)

### Common Score Improvements

| Refactoring Action                  | Typical Score Impact |
| ----------------------------------- | -------------------- |
| Split monolithic skill (>500 lines) | +10-15 points        |
| Standardize command frontmatter     | +5-8 points          |
| Resolve orphaned files              | +3-5 points          |
| Add missing plugin.json             | +5-7 points          |
| Fix broken reference links          | +5-10 points         |
| Optimize skill descriptions         | +3-5 points          |

### Quality Indicators

After refactoring, plugins should demonstrate:

- ✅ All skills <500 lines (preferably <400)
- ✅ No orphaned files
- ✅ All frontmatter validates against schema
- ✅ All reference links resolve
- ✅ Pre-commit hooks pass
- ✅ Skills load successfully via symlinks
- ✅ Commands appear in `/help` menu

---

## Git Integration

### Committing Refactoring Work

**Before Major Changes**:

```bash
git add .claude/plan/REFACTOR-PLAN.md .claude/plan/tasks-refactor-python3-development.md
git commit -m "plan: add python3-development refactoring tasks"
```

**During Refactoring Phases**:

```bash
# After Phase 1 (Foundation)
git add ./plugins/python3-development/
git commit -m "refactor(python3-dev): create skill directories and plugin.json"

# After Phase 2 (Skill Splits)
git add ./plugins/python3-development/skills/
git commit -m "refactor(python3-dev): split monolithic skill into 5 focused skills"

# After Phase 3 (Migration)
git add ./plugins/python3-development/skills/
git commit -m "refactor(python3-dev): migrate references and assets to python3-project"

# After Phase 4 (Validation)
git add ./plugins/python3-development/
git commit -m "refactor(python3-dev): complete refactoring (score 62→85)"
```

**Updating Plan Status**:

```bash
git add .claude/plan/REFACTOR-PLAN.md
git commit -m "plan: move python3-development to completed (score 62→85)"
```

---

## References

### Related Documentation

- [Plugin Assessment Criteria](../../plugins/holistic-linting/docs/assessment-criteria.md)
- [Skill Refactoring Guidelines](../../.claude/skills/skill-refactorer/SKILL.md)
- [Claude Skills Schema](../../.claude/skills/claude-skills-overview-2026/SKILL.md)
- [Plugin Structure Best Practices](../../.claude/skills/claude-plugins-reference-2026/SKILL.md)

### Related Agents

- `plugin-assessor` - Phase 1 assessment and Phase 5 validation
- `python-cli-design-spec` - Phase 2 refactoring design
- `task-planner` - Phase 3 task breakdown (this document)
- `skill-refactorer` - Phase 4 skill splitting execution
- `claude-context-optimizer` - Phase 4 frontmatter and doc optimization
- `plugin-docs-writer` - Phase 5 user documentation generation

---

## Document Metadata

- **Created**: 2026-01-22
- **Last Updated**: 2026-01-22
- **Version**: 1.0
- **Purpose**: Central index for tracking plugin refactoring projects across all phases
