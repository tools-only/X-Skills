# SAM Infrastructure Research and Design Prompt

## Mission

Research 5 existing agent coordination systems to understand their meta-messaging, artifact storage, and git worktree awareness patterns. Use findings to design a comprehensive infrastructure layer for the Stateless Agent Methodology (SAM).

## Background Context

**What is SAM?**: The Stateless Agent Methodology treats Claude as a stateless function operating in 7 stages:

1. Discovery & Interview
2. Context Gathering
3. Research & Learning
4. Design & Architecture
5. Planning
6. Implementation & Validation
7. Delivery

**Current Gap**: SAM currently lacks infrastructure for:

- Meta-messaging (agent-to-agent async communication via message queues)
- Task-agent association (YAML config for agent selection + skill loading per task/stage)
- Structured artifact storage (designs, specs, plans, tasks, checkpoints)
- MCP server interface to the infrastructure
- Git worktree awareness (pull, rebase, read state from commits)

**Reference Document**: Read `./methodology_development/stateless-agent-methodology.md` for complete SAM specification (relative to repo root).

## Target Repositories for Research

Research these 5 systems to understand their infrastructure patterns:

1. **~/repos/gastown** - Emphasize git worktree awareness patterns (user noted this has good example)
2. **~/repos/octocode-mcp** - MCP server integration patterns
3. **~/repos/get-shit-done** - Task and checkpoint management
4. **~/repos/cc-sessions** - Session state and continuity
5. **~/repos/BMAD-METHOD** - Multi-agent coordination (21+ agents)

## Research Focus Areas

For each repository, identify:

### 1. Meta-Messaging Patterns

**Definition**: Information about the process or tasks in the process - NOT code changes to the project.

**Questions to answer**:

- How do agents communicate asynchronously?
- What message queue systems are used (if any)?
- What types of meta-messages are sent between agents?
- How are messages structured (schema/format)?
- How do agents read from their message queues?
- Are messages persistent or ephemeral?

**Example meta-messages to look for**:

- "Task X completed, Task Y can now start"
- "Context gathering found 15 relevant files"
- "Design validation failed: missing constraint Y"
- "Agent A blocked on user clarification"

### 2. Artifact Types and Lifecycle

**Definition**: Structured documents that capture process state, decisions, plans, and progress.

**Questions to answer**:

- What artifact types are created? (designs, specs, plans, tasks, checkpoints, etc.)
- What is the lifecycle of each artifact? (created → updated → validated → archived)
- What information does each artifact contain?
- How are artifacts referenced between stages/agents?
- Are artifacts versioned? How?

**Example artifacts to identify**:

- Design documents (architecture decisions, constraints)
- Context documents (gathered information, research findings)
- Specification documents (requirements, acceptance criteria)
- Plan documents (task breakdown, dependencies)
- Task documents (individual work units with status)
- Checkpoint documents (progress snapshots, state recovery points)

### 3. Storage Patterns

**Questions to answer**:

- Where are artifacts stored? (filesystem, sqlite3, cloud db, git issues/milestones)
- What directory structure or schema is used?
- How are artifacts organized by project/session/stage?
- How are artifacts retrieved by agents?
- What indexing or search capabilities exist?

**Look for**:

- `.planning/` or similar directories
- Database schemas (sqlite/postgres)
- Git integration (issues, milestones, projects)
- Naming conventions for files/records

### 4. Agent Coordination Patterns

**Questions to answer**:

- How are agents selected for tasks?
- How is agent-task association configured? (hardcoded, YAML, dynamic)
- What skills/capabilities must agents load per task?
- How do agents signal completion or blocking?
- How are dependencies between agents handled?

**Look for**:

- Configuration files (YAML, JSON) mapping tasks to agents
- Agent capability declarations
- Task dependency specifications
- Completion signaling mechanisms

### 5. Git Worktree Awareness (Priority: Gastown)

**Questions to answer**:

- How does the system detect git state changes?
- When does pull/rebase happen in the workflow?
- How are commit messages used for state recovery?
- How does the system handle merge conflicts?
- How does the system read recent changes from git history?
- What git operations are first-class in the workflow?

**Look for** (especially in gastown):

- Git operation automation
- Commit message parsing for state
- Rebase handling
- Working tree status detection
- Change detection between sessions

## Output Requirements

Create two documents:

### Document 1: Infrastructure Research Findings

**Path**: `./methodology_development/sam-infrastructure-research.md` (relative to repo root)

**Structure**:

```markdown
# SAM Infrastructure Research Findings

## Executive Summary
- Key patterns identified across 5 systems
- Common approaches vs divergent approaches
- Recommended patterns for SAM

## Repository Analysis

### 1. Gastown
#### Meta-Messaging
[findings]

#### Artifact Storage
[findings]

#### Agent Coordination
[findings]

#### Git Worktree Awareness
[detailed findings - this is the priority example]

### 2. OctoCode MCP
[repeat structure]

### 3. Get Shit Done
[repeat structure]

### 4. CC-Sessions
[repeat structure]

### 5. BMAD-METHOD
[repeat structure]

## Cross-System Comparison

### Meta-Message Types (Table)
| System | Message Types | Queue Mechanism | Persistence |
|--------|---------------|-----------------|-------------|
| gastown | ... | ... | ... |

### Artifact Types (Table)
| System | Artifact Types | Storage | Versioning |
|--------|----------------|---------|------------|
| gastown | ... | ... | ... |

### Backend Storage (Table)
| System | Backend | Why Chosen | Tradeoffs |
|--------|---------|------------|-----------|
| gastown | ... | ... | ... |

## Pattern Recommendations for SAM

### Meta-Messaging Recommendation
[synthesis]

### Artifact Storage Recommendation
[synthesis]

### Agent Coordination Recommendation
[synthesis]

### Git Worktree Integration Recommendation
[synthesis]

## Evidence Log
[cite all sources with line numbers/file paths]
```

### Document 2: SAM Infrastructure Layer Design

**Path**: `./methodology_development/sam-infrastructure-layer.md` (relative to repo root)

**Structure**:

```markdown
# SAM Infrastructure Layer Specification

## Overview
High-level description of the infrastructure layer that connects SAM stages.

## Architecture

### Component 1: Meta-Messaging System
- Message queue design
- Message schema
- Agent subscription model
- Message persistence strategy

### Component 2: YAML Configuration System
- Task-agent association format
- Skill loading specification
- Stage-level configuration
- Example configurations

### Component 3: Artifact Storage System
- Artifact schema for each type (design, spec, plan, task, checkpoint)
- Storage backend options (filesystem, sqlite, supabase, git)
- Directory structure or database schema
- Retrieval API

### Component 4: MCP Server Interface
- Tool definitions for artifact CRUD
- Resource definitions for artifact access
- Prompt definitions for SAM stage guidance
- Integration with existing SAM skills/agents

### Component 5: Git Worktree Integration
- Auto-pull/rebase workflow
- Commit message parsing for state recovery
- Working tree status monitoring
- Conflict resolution protocol

## Implementation Roadmap
Phased approach to building the infrastructure layer.

## Backend Options Analysis
Comparison of filesystem vs sqlite vs cloud vs git-based backends.

## Integration with Existing SAM
How this infrastructure plugs into the 7-stage SAM pipeline.
```

## Research Methodology

Follow this process:

1. **For each repository**:

   - Start by reading README.md or main documentation
   - Use `Glob` to find configuration files (_.yaml,_.json, \*.toml)
   - Use `Glob` to find artifact storage directories (.planning/, .artifacts/, etc.)
   - Use `Grep` to search for message queue patterns, agent coordination keywords
   - Read key implementation files that handle messaging, artifacts, git operations
   - Document findings with file paths and line numbers for citations

2. **Evidence-based analysis**:

   - Cite specific files and line numbers for all claims
   - If a pattern is not explicitly found, mark as "NOT_OBSERVED" rather than assuming
   - Distinguish between documented patterns and implemented patterns
   - Note version/date of repository when accessed

3. **Synthesis**:
   - Create comparison tables showing patterns across all 5 systems
   - Identify common patterns (convergence) and unique approaches (divergence)
   - Recommend patterns for SAM based on evidence, not assumptions

## Tools Available

- `Read(file_path)` - Read any file in the repositories
- `Glob(pattern, path)` - Find files matching patterns
- `Grep(pattern, path, output_mode)` - Search file contents
- `Bash` - Run git commands or directory listing if needed

## Success Criteria

You have completed this task when:

- [ ] All 5 repositories have been researched with evidence citations
- [ ] Meta-message types documented for each system with examples
- [ ] Artifact types documented for each system with schemas/examples
- [ ] Storage patterns documented with directory structures or DB schemas
- [ ] Git worktree awareness patterns documented (especially gastown)
- [ ] Comparison tables created showing patterns across systems
- [ ] sam-infrastructure-research.md created with complete findings
- [ ] sam-infrastructure-layer.md created with design specification
- [ ] All claims cited with file paths and line numbers
- [ ] Recommendations grounded in observed patterns, not assumptions

## Starting Point

Begin with gastown repository since user emphasized it has good git worktree awareness example:

```
cd ~/repos/gastown
Read README.md or similar entry point
Glob **/*.{yaml,json,toml} to find config files
Grep for git-related keywords (rebase, commit, worktree, etc.)
```

## Context Files to Read First

Before starting research, read these files to understand SAM:

1. `./methodology_development/stateless-agent-methodology.md` - Complete SAM specification
2. `./methodology_development/sam-framework-generator.md` - Meta-SAM concept (using SAM to generate SAM implementations)

## Notes

- Focus on **meta-messages** (information about process/tasks), not code change messages
- Git worktree awareness is critical - gastown is the priority example
- Evidence-based research only - cite all findings with file paths and line numbers
- If uncertain, mark "UNVERIFIED" and explain the gap rather than assuming
