# Backlog: Plugin-Creator Refactor Workflow Disconnects

## Summary

The plugin-creator refactor workflow (assessor → implement-refactor → start-refactor-task → ensure-complete) has 4 verified disconnects between steps. Identified by reading source files on 2026-02-09.

## Disconnect Inventory

### D1: Assessor TaskCreate entries invisible to sub-agents

**Files**: `plugins/plugin-creator/skills/assessor/SKILL.md`
**Evidence**: Lines 26-44 create 6 TaskCreate entries. Sub-agent prompts at lines 51-134, 175-336, 377-552, 600-646 do not pass task IDs.
**Observed**: Tasks #7-12 created during holistic-linting refactor were never updated past Phase 1 despite successful completion.

### D2: Inconsistent tracking API across workflow

**Files**:
- `plugins/plugin-creator/skills/assessor/SKILL.md` — uses TaskCreate (lines 29-34)
- `plugins/plugin-creator/skills/implement-refactor/SKILL.md` — uses TodoWrite (lines 54-59)
- `plugins/plugin-creator/skills/ensure-complete/SKILL.md` — uses TodoWrite (lines 27-33)

**Observed**: Two different APIs for the same purpose within one workflow.

### D3: Missing completion protocol instruction

**Files**: `plugins/plugin-creator/skills/implement-refactor/SKILL.md`
**Evidence**: Lines 109-116 check for task status ✅ COMPLETE in persistent task file, but do not instruct sub-agents to call `/start-refactor-task --complete {id}`.
**Related**: `plugins/plugin-creator/skills/start-refactor-task/SKILL.md` lines 26-31 define the `--complete` mechanism, line 150 shows self-invocation pattern.

### D4: Session-scoped TaskCreate entries don't survive context compaction

**Files**: `plugins/plugin-creator/skills/assessor/SKILL.md`
**Evidence**: Lines 29-34 use TaskCreate (session memory only). No Write tool calls persist these to files. Long-running assessor sessions lose progress tracking when context is compacted.
**Observed**: Holistic-linting refactor completed (score 68→92) but TaskCreate entries showed stale state.

## Tracking Systems in Use

| System | Scope | Created By | Persists? |
|--------|-------|------------|-----------|
| TaskCreate/TaskUpdate | Session memory | assessor, start-refactor-task | No |
| TodoWrite | Session memory | implement-refactor, ensure-complete | No |
| `.claude/plan/tasks-refactor-{slug}.md` | File | swarm-task-planner (via assessor Phase 3) | Yes |
| `.claude/plan/REFACTOR-PLAN.md` | File | swarm-task-planner, ensure-complete | Yes |

## Source Investigation

- Agent run on 2026-02-09, read 5 skill files (2,074 lines), 2 plan files (520 lines)
- Skills read: refactor-plugin, assessor, implement-refactor, ensure-complete, start-refactor-task
- All line references verified against source files

## Status

**Not started** — created as backlog item for separate session.
