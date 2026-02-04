# SAM Harness

**Date**: 2026-01-27
**Status**: Working Document
**Purpose**: Define the implementation architecture for the Stateless Agent Methodology (SAM) as a Claude Code Plugin

---

## Overview

The **SAM Harness** is the AI workflow automation system that implements the Stateless Agent Methodology (SAM) as a Claude Code Plugin. It provides the orchestration, artifact management, and execution infrastructure to operationalize SAM's constraint-driven development framework.

## Architecture

```
Stateless Agent Methodology (SAM)
    ↓ (methodology defines)
SAM Harness
    ↓ (packaged as)
Claude Code Plugin
    ├── Scripts (orchestration, task management, git worktrees)
    ├── Skills (SAM stage workflows)
    ├── Agents (specialized workers for each stage)
    └── MCP Server (artifact storage, state management)
```

## Components

### 1. Scripts

**Purpose**: Orchestration logic and infrastructure management

**Responsibilities**:

- Pipeline orchestration (stage transitions)
- Task queue management
- Git worktree creation and cleanup
- Parallel/sequential execution coordination
- Progress tracking and status reporting

**Example scripts**:

- `orchestrator.py` - Main pipeline coordinator
- `worktree-manager.py` - Git worktree lifecycle
- `task-dispatcher.py` - Worker assignment and execution

### 2. Skills

**Purpose**: SAM stage workflow definitions

**Mapping to SAM Stages**:

| Skill                 | SAM Stage                    | Input               | Output                                      |
| --------------------- | ---------------------------- | ------------------- | ------------------------------------------- |
| `discovery`           | Stage 1: Discovery           | User request        | `ARTIFACT:DISCOVERY(SCOPE:...)`             |
| `planning`            | Stage 2: Planning (RT-ICA)   | Discovery artifacts | `ARTIFACT:PLAN(SCOPE:...)`                  |
| `context-integration` | Stage 3: Context Integration | Plan + codebase     | `ARTIFACT:PLAN(SCOPE:...)` (contextualized) |
| `task-decomposition`  | Stage 4: Task Decomposition  | Contextualized plan | `ARTIFACT:TASK(TASK:...)` files             |
| `execution`           | Stage 5: Execution           | Single task file    | `ARTIFACT:EXECUTION(TASK:...)`              |
| `forensic-review`     | Stage 6: Forensic Review     | Execution results   | `ARTIFACT:REVIEW(TASK:...)`                 |
| `final-verification`  | Stage 7: Final Verification  | All completed tasks | `ARTIFACT:VERIFICATION(SCOPE:...)`          |

### 3. Agents

**Purpose**: Specialized workers invoked by skills for parallel/focused work

**Agent Types**:

| Agent                  | Invoked By                     | Specialization                                   |
| ---------------------- | ------------------------------ | ------------------------------------------------ |
| `research-agent`       | Discovery skill                | Online research, package documentation lookup    |
| `codebase-explorer`    | Discovery, Context Integration | Local repository analysis, git forensics         |
| `prerequisite-checker` | Planning skill                 | RT-ICA assessment, dependency verification       |
| `task-creator`         | Task Decomposition skill       | Atomic task generation with embedded context     |
| `implementation-agent` | Execution skill                | Code implementation following task specification |
| `quality-checker`      | Forensic Review skill          | Independent verification, fact-checking          |

**Key Constraint**: Each agent executes exactly one task, then terminates (no multi-task agents).

### 4. MCP Server

**Purpose**: Artifact storage and state management

**Responsibilities**:

- CRUD operations for `ARTIFACT:{TYPE}(SCOPE:...)` tokens
- Artifact versioning and history
- Task queue state persistence
- Execution status tracking
- Query interface for artifact retrieval

**Storage Backend Options**:

- Filesystem-backed (`.sam/artifacts/`)
- SQLite database
- PostgreSQL (for multi-user/concurrent workflows)

**API Endpoints** (example):

- `create_artifact(token, content)` → Store new artifact
- `read_artifact(token)` → Retrieve artifact content
- `update_artifact(token, content)` → Update existing artifact
- `list_artifacts(type, scope)` → Query artifacts by type/scope
- `get_task_status(task_id)` → Retrieve task execution state

## Plugin Structure

```
.claude/plugins/sam-harness/
├── plugin.json                 # Plugin metadata and configuration
├── scripts/
│   ├── orchestrator.py         # Pipeline coordinator
│   ├── worktree-manager.py     # Git worktree lifecycle
│   └── task-dispatcher.py      # Worker assignment
├── skills/
│   ├── discovery/
│   │   └── SKILL.md
│   ├── planning/
│   │   └── SKILL.md
│   ├── context-integration/
│   │   └── SKILL.md
│   ├── task-decomposition/
│   │   └── SKILL.md
│   ├── execution/
│   │   └── SKILL.md
│   ├── forensic-review/
│   │   └── SKILL.md
│   └── final-verification/
│       └── SKILL.md
├── agents/
│   ├── research-agent.md
│   ├── codebase-explorer.md
│   ├── prerequisite-checker.md
│   ├── task-creator.md
│   ├── implementation-agent.md
│   └── quality-checker.md
└── mcp-server/
    ├── server.py               # MCP server implementation
    ├── storage/
    │   ├── filesystem.py       # Filesystem backend
    │   └── sql.py              # SQL backend
    └── schema.sql              # Database schema (if using SQL)
```

## Workflow Execution

### User Initiates Request

```bash
# User invokes SAM workflow
/sam start "Add health check endpoint"
```

### Harness Execution Flow

1. **Orchestrator** receives request
2. **Discovery skill** activates, dispatches research agents
3. Discovery produces `ARTIFACT:DISCOVERY(SCOPE:health-check)`
4. **Planning skill** activates, runs RT-ICA gate
5. Planning produces `ARTIFACT:PLAN(SCOPE:health-check)`
6. **Context Integration skill** activates, explores codebase
7. Context produces `ARTIFACT:PLAN(SCOPE:health-check)` (contextualized)
8. **Task Decomposition skill** activates, creates atomic tasks
9. Task Decomposition produces `ARTIFACT:TASK(TASK:001)`, `ARTIFACT:TASK(TASK:002)`, etc.
10. **Orchestrator** dispatches execution agents (parallel/sequential based on dependencies)
11. Each **Execution agent** (fresh session) executes one task
12. Execution produces `ARTIFACT:EXECUTION(TASK:001)`, etc.
13. **Forensic Review agents** verify each execution
14. Review produces `ARTIFACT:REVIEW(TASK:001)`, etc.
15. **Orchestrator** loops if NEEDS_WORK found
16. **Final Verification skill** certifies feature completion
17. Final Verification produces `ARTIFACT:VERIFICATION(SCOPE:health-check)`

### State Persistence

All artifacts stored via MCP server:

- Survives session resets
- Queryable for progress tracking
- Enables resume-from-checkpoint
- Audit trail for debugging

## Key Design Decisions

### 1. Scripts vs Skills vs Agents

| Component   | Purpose                       | Language | Invocation                          |
| ----------- | ----------------------------- | -------- | ----------------------------------- |
| **Scripts** | Infrastructure, orchestration | Python   | Direct execution by orchestrator    |
| **Skills**  | Workflow stage definitions    | Markdown | Activated by user or orchestrator   |
| **Agents**  | Specialized workers           | Markdown | Invoked via `Task()` tool by skills |

### 2. MCP Server as Single Source of Truth

- All state externalized to MCP server
- No state in conversation history
- No state in agent memory
- Artifacts are durable, versioned, queryable

### 3. One Task Per Agent

- Execution agents are single-use
- Fresh session per task (no context accumulation)
- Eliminates error propagation between tasks
- Enforces statelessness principle

### 4. Git Worktrees for Parallel Execution

- Each task executes in isolated worktree
- Prevents file conflicts during parallel work
- Clean merging after verification passes
- Automatic cleanup on completion

## Integration Points

### With SAM Methodology

The harness implements SAM's 7-stage pipeline:

- Enforces RT-ICA gate (blocks on missing prerequisites)
- Ensures task files contain complete context
- Runs independent forensic review (not self-verification)
- Implements artifact-based state transitions

### With Claude Code

The harness extends Claude Code via:

- Plugin system (installable, versioned)
- Skills (user-invocable workflows)
- Agents (sub-agent delegation via `Task()`)
- MCP server (tool integration for artifact management)

### With Existing Development Tools

The harness integrates with:

- Git (worktrees, commits, branches)
- Linters (deterministic backpressure)
- Test runners (verification gates)
- CI/CD (final deployment pipeline)

## Open Questions

1. **Naming**: Is "SAM Harness" the final name, or placeholder?
2. **MCP Server**: Filesystem or SQL backend by default?
3. **Orchestrator**: Pure Python script or hybrid with Claude orchestration?
4. **Task Representation**: How are `ARTIFACT:TASK(TASK:...)` files structured?
5. **Parallel Execution**: Max concurrent workers? Resource limits?
6. **Error Handling**: What happens when RT-ICA blocks? User intervention protocol?
7. **Resume Capability**: Can users pause/resume workflows? Checkpoint granularity?

## Next Steps

- [ ] Define plugin.json schema
- [ ] Implement MCP server prototype (filesystem backend)
- [ ] Create discovery skill with research agent delegation
- [ ] Build orchestrator script for stage transitions
- [ ] Test end-to-end workflow with simple feature request
- [ ] Document user-facing commands (`/sam start`, `/sam status`, etc.)
- [ ] Create installation and configuration guide

---

**Related Documentation**:

- [Stateless Agent Methodology](./stateless-agent-methodology.md)
- [Stateless Software Engineering Framework](./stateless-software-engineering-framework.md)
- [Process Realignment](./process_realignment.md)
