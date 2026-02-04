# Stateless Agent Methodology vs cc-sessions

A comparison of two approaches to improving AI pair programming with Claude Code.

---

## Overview

| Aspect             | Stateless Agent Methodology (SAM)                                                                        | cc-sessions                                                                                    |
| ------------------ | -------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------- |
| **Type**           | Methodology/process framework                                                                            | Workflow enforcement framework                                                                 |
| **Focus**          | Addressing LLM cognitive limitations through architecture                                                | Enforcing structured approval workflow before code changes                                     |
| **Implementation** | Conceptual framework requiring custom implementation                                                     | Ready-to-use npm/pip package with hooks and protocols                                          |
| **Target**         | Any LLM agent system                                                                                     | Claude Code specifically                                                                       |
| **Core mechanism** | Phase separation with artifact passing                                                                   | DAIC (Discussion-Alignment-Implementation-Check) enforcement                                   |
| **Backpressure**   | Deterministic checks as ground truth (tests/lint/static analysis/checklists); iterate until pass/BLOCKED | Tool gating + protocol enforcement; can incorporate checks but not the primary control surface |

---

## Core Philosophy

### Stateless Agent Methodology

**Central insight**: Claude is not a knowledge worker—Claude is a stateless computation engine. LLM agents cannot reliably self-assess knowledge gaps. They optimize for _apparent completion_ over _correct completion_.

**Solution**: Externalize state and enforcement into artifacts + gates:

- Separate stages with different agents
- Task files as the _complete prompt_ for execution
- Deterministic backpressure (tests/lint/static analysis/checklists) treated as ground truth
- Independent forensic review (not self-review)

**Key principle**: Treat Claude like a pure function: input complete context (task file with all answers), output verified result, no memory between invocations.

### cc-sessions

**Central insight**: Claude will start writing code immediately without discussion, leading to unwanted changes and scope creep.

**Solution**: Block write-capable tools by default. Claude must discuss, get explicit approval via trigger phrases, then implement only the approved todos.

**Key principle**: Claude earns the right to write code through explicit user approval of specific plans.

---

## Alignment on Core Problems

Both frameworks identify similar issues but address them differently:

| Problem                                        | SAM Approach                                                                   | cc-sessions Approach                                                                        |
| ---------------------------------------------- | ------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------- |
| **Premature implementation**                   | Assessment phase blocks until prerequisites verified                           | DAIC blocks tools until explicit approval                                                   |
| **Scope creep**                                | Task files contain exact scope; forensic review catches deviations             | Todo validation locks approved plan; changes trigger "SHAME RITUAL"                         |
| **Long-context degradation (“context rot”)**   | Avoid long sessions; keep contexts bounded via task files and stage separation | Manages long sessions via protocols (discussion/approval) and context management/compaction |
| **Training data staleness (knowledge cutoff)** | Assume priors are stale; force grounding + verification via artifacts/tools    | Encourages investigation-first workflow; can enforce process but not replace grounding      |
| **Context loss across sessions**               | Execution does not require conversation history; artifacts are source of truth | Persistent task/runtime state enables resuming the same workflow                            |
| **Self-verification bias**                     | Separate forensic review agent (independent)                                   | Same session validates, but via prescribed formats and todo tracking                        |
| **Methodology skipping**                       | Methodology IS the prompt structure                                            | Hook-based enforcement blocks tools at system level                                         |

---

## Architecture Comparison

### Enforcement Mechanism

**SAM**: Structural enforcement through pipeline architecture

```
User Request → Discovery → Planning (RT-ICA gate) → Context Integration
→ Task Decomposition → Execution → Forensic Review → (Orchestration Loop) → Final Verification
```

**cc-sessions**: Hook-based enforcement through DAIC workflow

```
Discussion Mode (tools blocked) → User approval trigger →
Implementation Mode (tools enabled) → Todo completion →
Return to Discussion Mode
```

### Phase/Protocol Mapping

| SAM Stage                        | cc-sessions Equivalent                        |
| -------------------------------- | --------------------------------------------- |
| **Stage 1: Discovery**           | Discussion mode + task-creation protocol      |
| **Stage 2: Planning (RT-ICA)**   | DAIC alignment phase + todo proposal          |
| **Stage 3: Context Integration** | context-gathering agent                       |
| **Stage 4: Task Decomposition**  | Task creation with markdown frontmatter       |
| **Stage 5: Execution**           | Implementation mode with todo tracking        |
| **Stage 6: Forensic Review**     | code-review agent (but in same orchestration) |
| **Stage 7: Final Verification**  | task-completion protocol                      |
| **(Orchestration Loop)**         | Automatic mode transitions                    |

### Agent Separation

**SAM** separates agents by _function_:

- Discovery agent (questioning, gathering)
- Planning agent (RT-ICA, design)
- Context Integration agent (codebase mapping)
- Task Decomposition agent (atomic task creation)
- Execution agent (implementation)
- Forensic Review agent (independent verification)
- Final Verification agent (goal validation)

**cc-sessions** separates agents by _operation type_:

- context-gathering (creates task context manifests)
- logging (consolidates work logs)
- code-review (reviews implementations)
- context-refinement (updates context with discoveries)
- service-documentation (maintains CLAUDE.md files)

**Key difference**: SAM agents are pipeline stages. cc-sessions agents are utility workers called from main session.

---

## Key Differences

### 1. Statelessness vs State Persistence

| Aspect                   | SAM                                             | cc-sessions                                 |
| ------------------------ | ----------------------------------------------- | ------------------------------------------- |
| **Philosophy**           | Intentionally stateless—fresh context per agent | Persistent state across sessions            |
| **Session continuity**   | Task file IS the complete context               | sessions-state.json preserves runtime state |
| **Cross-session memory** | None by design                                  | Task state, todos, mode all persist         |
| **Restart behavior**     | Agent receives task file anew                   | Resumes from saved state                    |

**SAM**: Each agent is a pure function. Task file contains everything needed. No conversation history.

**SAM**: Stateless _sessions_, persistent _artifacts_. Each execution agent is a pure function over a task file + repo state; durability comes from writing/reading artifacts, not from conversation history.

**cc-sessions**: State persists in `sessions-state.json`. Resume where you left off. Mode, todos, active protocol all survive restarts.

### 2. Verification Approach

**SAM**: Independent forensic review agent

- Different agent than executor
- Cannot see execution conversation
- Validates against task specification
- Result: COMPLETE or NEEDS_WORK

**cc-sessions**: Same session validation

- code-review agent for implementations
- Todo completion detection
- Prescribed output formats ([PROPOSAL], [STATUS], etc.)
- "SHAME RITUAL" for unauthorized changes

### 3. Approval Mechanism

**SAM**: RT-ICA gate (prerequisite-based)

- Lists all prerequisites
- Marks each: AVAILABLE | DERIVABLE | MISSING
- BLOCKS if any MISSING
- Structural gate—cannot proceed without prerequisites

**cc-sessions**: DAIC (approval-based)

- Discussion mode default (tools blocked)
- Trigger phrases unlock implementation ("yert", "make it so")
- Todo list defines exact scope
- Changes to approved todos blocked with diff display

### 4. Tool Blocking

**SAM**: No explicit tool blocking—agents receive only appropriate tools for their phase

**cc-sessions**: Hook-based tool blocking

- Edit, Write, MultiEdit blocked in discussion mode
- 70+ read-only bash commands recognized
- Write operations detected in complex pipelines
- Configurable blocked tool lists

### 5. Context Management

**SAM**: Complete context in task file

- All answers provided
- No training data recall needed
- File/URL references explicit
- Definition of Done embedded

**cc-sessions**: Context gathering via agent

- context-gathering agent creates manifests
- context-refinement updates with discoveries
- Context compaction protocol for long sessions
- "squish" trigger for manual compaction

---

## Feature Comparison Matrix

| Feature                   | SAM                                  | cc-sessions                                                                       |
| ------------------------- | ------------------------------------ | --------------------------------------------------------------------------------- |
| Tool blocking             | N/A (structural separation)          | Hook-based DAIC enforcement                                                       |
| Prerequisite gates        | RT-ICA (AVAILABLE/DERIVABLE/MISSING) | Discussion mode approval                                                          |
| Scope control             | Task file defines exact scope        | Todo list locked after approval                                                   |
| Change detection          | Forensic review agent                | Todo validation with diff display                                                 |
| Cross-session persistence | None (stateless by design)           | Full state in sessions-state.json                                                 |
| Git integration           | Not specified                        | Git-aware workflow enforcement/protocols (verify exact features per installation) |
| Trigger phrases           | N/A                                  | Customizable ("yert", "mek:", etc.)                                               |
| Specialized agents        | 7 pipeline stages                    | 5 utility agents                                                                  |
| Protocol automation       | Artifact flow between stages         | Protocol templates with variables                                                 |
| Installation              | Conceptual (implement yourself)      | npm/pip package with installer                                                    |

---

## Complementary Strengths

These frameworks address overlapping but distinct concerns:

| Layer            | SAM                                               | cc-sessions                              |
| ---------------- | ------------------------------------------------- | ---------------------------------------- |
| **Cognitive**    | Phase separation prevents self-assessment failure | DAIC prevents premature implementation   |
| **Enforcement**  | Structural (methodology IS the prompt)            | Behavioral (hook-based tool blocking)    |
| **Verification** | External forensic agent                           | Same-session code-review + todo tracking |
| **Persistence**  | None (intentional)                                | Full state preservation                  |
| **Tooling**      | Conceptual framework                              | Production-ready package                 |

### Where cc-sessions Implements SAM Patterns

cc-sessions already implements several SAM concepts:

1. **Prerequisite gates**: Discussion mode blocks until approval
2. **Scope control**: Todo list defines exact implementation scope
3. **Change detection**: Unauthorized todo changes blocked
4. **Specialized agents**: Heavy operations delegated to separate contexts
5. **Context manifests**: context-gathering creates task context files

### Where SAM Differs

1. **True statelessness**: Each agent is a fresh session
2. **External verification**: Forensic agent cannot see execution
3. **Complete context injection**: No session memory needed
4. **Artifact passing**: Communication via files, not shared context

---

## Integration Opportunities

### cc-sessions Could Adopt

1. **External forensic review**: Separate verification agent that doesn't see execution
2. **RT-ICA gate**: Block planning until prerequisites formally verified
3. **Stateless execution**: Option to spawn fresh agent for each task
4. **Artifact-based handoffs**: Reduce dependency on session state

### SAM Could Adopt

1. **Hook-based enforcement**: System-level tool blocking
2. **Trigger phrase system**: Natural language protocol activation
3. **Git integration**: Branch enforcement, commit automation
4. **State persistence**: Optional state preservation for long projects

---

## When to Use Each

### Use Stateless Agent Methodology When

- Maximum assurance required—no tolerance for hallucination
- Independent verification is critical
- Working with internal/proprietary codebases
- Context window pressure is a major concern
- Building custom agent orchestration systems
- Need provable verification trail

### Use cc-sessions When

- Need ready-to-use enforcement immediately
- Working in Claude Code environment
- Want approval-based workflow (discuss before implement)
- Need task persistence across session restarts
- Want git workflow automation (branches, commits, merges)
- Prefer trigger phrase control ("yert", "finito")

### Use Both When

- Want cc-sessions' enforcement with SAM's external verification
- Need persistent tasks but stateless execution
- Building high-assurance features with practical tooling
- Implementing phased workflow on top of cc-sessions protocols

---

## Summary

**Stateless Agent Methodology** is a _cognitive architecture_ that treats Claude as a pure function—input complete context, output verified result, no memory between invocations. Verification is external and independent.

**cc-sessions** is a _workflow enforcement system_ that requires Claude to earn the right to write code through explicit approval. It persists state across sessions and uses hook-based tool blocking to enforce the Discussion-Alignment-Implementation-Check pattern.

**Key philosophical differences**:

| Dimension            | SAM                           | cc-sessions                    |
| -------------------- | ----------------------------- | ------------------------------ |
| **Trust model**      | Don't trust—verify externally | Trust after explicit approval  |
| **State philosophy** | Stateless is a feature        | Persistence enables continuity |
| **Verification**     | Different agent than executor | Same session, different phase  |
| **Enforcement**      | Structural (pipeline)         | Behavioral (hooks)             |

**Both agree on**:

- Claude will skip methodology if not enforced
- Explicit approval before implementation is essential
- Scope creep must be actively prevented
- Structured workflows improve reliability
- Context management is critical

**The most robust approach** may be implementing SAM's external forensic review within cc-sessions' infrastructure—combining the practical tooling and state persistence of cc-sessions with the independent verification of SAM.

---

## References

- [Stateless Agent Methodology](./stateless-agent-methodology.md) - Full methodology documentation
- [cc-sessions GitHub](https://github.com/GWUDCAP/cc-sessions) - Official repository
- [cc-sessions CLAUDE.md](https://github.com/GWUDCAP/cc-sessions/blob/main/CLAUDE.md) - Architecture documentation
- [cc-sessions CHANGELOG](https://github.com/GWUDCAP/cc-sessions/blob/main/CHANGELOG.md) - Feature history
