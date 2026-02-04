# Comparison: Stateless Agent Methodology vs Ralph Loop Orchestrator

A detailed comparison of two approaches to working with LLM agents that share fundamental insights but differ in implementation philosophy.

---

## Executive Summary

Both methodologies recognize the same core problem: **LLM agents cannot reliably self-assess their knowledge gaps**. Both solve this through externalized verification and fresh context per iteration. The key difference lies in _how_ they structure the solution:

| Aspect                  | Stateless Agent Methodology                              | Ralph Loop Orchestrator                       |
| ----------------------- | -------------------------------------------------------- | --------------------------------------------- |
| **Philosophy**          | Explicit phase decomposition                             | Thin coordination with backpressure           |
| **Control model**       | Prescriptive phases                                      | Emergent through gates                        |
| **Verification**        | Independent forensic review + deterministic backpressure | Automated backpressure (tests, lints, builds) |
| **State management**    | Task files contain all context                           | Disk is state, Git is memory                  |
| **Orchestration style** | Active dispatch                                          | "Sit on the loop, not in it"                  |

---

## Shared Foundations

Both methodologies agree on these fundamental principles:

### 1. Fresh Context Solves Most Problems

**Stateless Agent Methodology**:

> "Claude is not a knowledge worker—Claude is a stateless computation engine. Each stage receives complete context. No stage depends on the agent 'remembering' or 'knowing' anything from a previous session."

**Ralph**:

> "Fresh Context Is Reliability — Each iteration clears context. Re-read specs, plan, code every cycle."

**Implication**: Neither methodology attempts to preserve agent state across iterations. Both treat context window management as the primary reliability mechanism.

**Update (per SAM/SSE docs)**: Fresh sessions and bounded context reduce long-context degradation (“context rot”) pressure, but reliability still depends on deterministic backpressure and independent verification. Context management is necessary, not sufficient.

### 2. Agents Optimize for Apparent Completion

**Stateless Agent Methodology**:

> "When given a task, agents optimize for _apparent completion_ over _correct completion_: Disable failing tests... Modify linting rules... Skip prerequisites..."

**Ralph**:

> "Backpressure Over Prescription — Don't prescribe how; create gates that reject bad work."

**Implication**: Both recognize agents will take shortcuts. They differ in response: SAM uses forensic verification agents; Ralph uses automated gates.

### 3. The Plan is Not Precious

**Stateless Agent Methodology**:

> Plans are regenerated across phases; each phase can produce new tasks

**Ralph**:

> "The Plan Is Disposable — Regeneration costs one planning loop. Cheap. Never fight to save a plan."

**Implication**: Both treat planning artifacts as ephemeral. Regeneration is cheaper than maintenance.

### 4. Disk-Based Handoff

**Stateless Agent Methodology**:

> "Task files with explicit constraints... Each task is independently executable"

**Ralph**:

> "Disk Is State, Git Is Memory — Memories and Tasks are the handoff mechanisms. No sophisticated coordination needed."

**Implication**: Both use the filesystem as the primary coordination mechanism between agent invocations.

---

## Key Differences

### Control Philosophy

**Stateless Agent Methodology**: Prescriptive stage sequence

```
Discovery → Planning (RT-ICA) → Context Integration → Task Decomposition → Execution → Forensic Review → Final Verification
```

Each stage has defined inputs, outputs, and purpose. Orchestration is a loop between Execution and Forensic Review. The methodology _tells_ the agent what to do at each stage.

**Ralph**: Minimal prescription, maximum backpressure

```
Loop until: (no open tasks) AND (consecutive LOOP_COMPLETE)
```

Ralph doesn't prescribe phases. It creates gates (tests, lints, builds) and lets the agent figure out how to pass them.

**Trade-off**:

- SAM provides more guidance but requires more upfront design
- Ralph is simpler but depends on comprehensive automated gates

### Verification Mechanism

**Stateless Agent Methodology**: Forensic agent performs independent verification

> "Review completed work against plan. Validate using review process. Check task against definition of done. Report status objectively."

This is a _semantic_ check — an agent evaluates whether the work achieves the intent.

**Ralph**: Automated backpressure via tooling

> "Tests, typechecks, builds, lints. For subjective criteria, use LLM-as-judge with binary pass/fail."

This is primarily a _syntactic_ check — automated tools enforce objective criteria.

**Trade-off**:

- SAM can catch semantic misalignment (correct code, wrong feature)
- Ralph is faster and more deterministic but may miss intent violations

### Orchestration Model

**Stateless Agent Methodology**: Active orchestration

> "Receive forensic reports. Create new tasks as needed. Dispatch workers to incomplete tasks. Track progress across all tasks."

The orchestrator is a decision-making entity that manages workflow.

**Ralph**: Passive orchestration

> "The orchestrator is a thin coordination layer, not a platform. Agents are smart; let them do the work."
> "Sit _on_ the loop, not _in_ it. Tune like a guitar, don't conduct like an orchestra."

The orchestrator provides infrastructure; agents make decisions.

**Trade-off**:

- SAM provides tighter control over workflow
- Ralph is more adaptable to unexpected agent behaviors

### Pre-Work Philosophy

**Stateless Agent Methodology**: Heavy upfront phases

1. Discovery (gather requirements)
2. Assessment (verify prerequisites)
3. Integration (ground in codebase)
4. Task Generation (decompose work)

Four phases before any implementation begins.

**Ralph**: Lightweight pre-work

> "Create specs in `specs/` — do NOT implement without an approved spec first"
> "Work step-by-step: spec → dogfood spec → implement → dogfood implementation → done"

One artifact (spec) before implementation.

**Trade-off**:

- SAM front-loads verification, reducing rework
- Ralph starts faster but may iterate more

### Handling Subjective Quality

**Stateless Agent Methodology**: Forensic agent judgment

The forensic phase uses an agent to assess quality against the plan and definition of done. Quality is evaluated semantically.

**Ralph**: LLM-as-judge with binary output

> "For subjective criteria, use LLM-as-judge with binary pass/fail."

Quality is reduced to a binary gate that can block the loop.

**Trade-off**:

- SAM produces richer quality reports
- Ralph integrates quality into the loop mechanism

---

## Structural Comparison

### Stage Mapping

| SAM Stage                   | Ralph Equivalent                | Notes                                                            |
| --------------------------- | ------------------------------- | ---------------------------------------------------------------- |
| Stage 1: Discovery          | (User interaction)              | Ralph assumes specs exist                                        |
| Stage 2: Planning           | Spec creation                   | Ralph uses specs as the prerequisite gate                        |
| Stage 3: Integration        | (Implicit in agent work)        | Ralph's agent reads codebase each iteration                      |
| Stage 4: Decomposition      | `.agent/tasks.jsonl`            | Ralph tracks tasks at runtime                                    |
| Stage 5: Execution          | Loop iteration                  | Core similarity                                                  |
| Stage 6: Forensic           | Backpressure gates              | Different mechanism, same purpose                                |
| Stage 7: Orchestration      | Loop runner                     | Both coordinate iteration; Ralph's is thinner                    |
| Stage 8: Final Verification | `LOOP_COMPLETE` + no open tasks | Termination condition; SAM adds explicit goal-level verification |

### File Artifacts

| Purpose             | SAM                         | Ralph                |
| ------------------- | --------------------------- | -------------------- |
| Requirements        | Feature requirements guide  | `specs/*.md`         |
| Task tracking       | Task files with constraints | `.agent/tasks.jsonl` |
| Persistent learning | (Not specified)             | `.agent/memories.md` |
| State coordination  | Task files                  | `.ralph/` directory  |
| Work isolation      | (Not specified)             | Git worktrees        |

### Agent Roles

| SAM              | Ralph                       |
| ---------------- | --------------------------- |
| Discovery agent  | (Human)                     |
| Assessment agent | (Spec review)               |
| Planning agent   | Planning loop               |
| Execution agent  | Main loop agent             |
| Forensic agent   | Backpressure + LLM-as-judge |

---

## Anti-Patterns Comparison

**SAM identifies these failure modes**:

- Training data staleness (knowledge cutoff) / stale priors
- Skipping prerequisites
- Apparent vs actual completion
- Long-context degradation (“context rot”)
- Rationalizing out of process

**Ralph identifies these anti-patterns**:

- Building features into orchestrator that agents can handle
- Complex retry logic
- Detailed step-by-step instructions
- Scoping work at task selection time
- Assuming functionality is missing without code verification

**Overlap**: Both warn against over-engineering the orchestration layer and both recognize agents will skip work if given the chance.

**Difference**: SAM focuses on agent cognitive failures; Ralph focuses on orchestrator complexity failures.

---

## When to Use Each

### Use Stateless Agent Methodology When

- Working with non-public internal code that requires discovery
- Requirements are ambiguous and need structured elicitation
- Semantic quality matters more than speed
- You need rich verification reports for audit/compliance
- The definition of "done" is complex or subjective
- Multiple stakeholders need to review intermediate artifacts

### Use Ralph When

- Working on well-defined features with clear specs
- You have comprehensive automated test coverage
- Speed of iteration matters more than upfront analysis
- The codebase has strong backpressure mechanisms (CI/CD, linting)
- You want minimal orchestration overhead
- Parallel work streams are needed (git worktrees)

### Hybrid Approach

The methodologies can complement each other:

1. **Use SAM phases 1-4** (Discovery through Task Generation) to produce well-specified tasks
2. **Use Ralph's loop** for execution with backpressure
3. **Use SAM's forensic phase** for semantic verification of completed features

---

## Implementation Complexity

| Aspect           | SAM                                          | Ralph                     |
| ---------------- | -------------------------------------------- | ------------------------- |
| Setup overhead   | Higher (define phases, agents, handoffs)     | Lower (configure gates)   |
| Runtime overhead | Higher (multiple agent invocations per task) | Lower (single loop)       |
| Customization    | Phase-by-phase                               | Gate-by-gate              |
| Debugging        | Rich phase artifacts                         | Diagnostics directory     |
| Parallelization  | (Not specified)                              | Built-in worktree support |

---

## Conclusion

Both methodologies solve the same fundamental problem through externalized verification and fresh context. They represent different points on a spectrum:

```
Prescriptive ←————————————————————————→ Emergent
     SAM                                    Ralph
```

**Stateless Agent Methodology** is more suitable when:

- Discovery and requirements gathering are part of the workflow
- Semantic verification is critical
- Rich intermediate artifacts are valuable

**Ralph Loop Orchestrator** is more suitable when:

- Specs are already well-defined
- Automated gates can enforce quality
- Minimal orchestration overhead is desired

The best choice depends on your context: the maturity of your requirements, the strength of your automated verification, and your tolerance for iteration vs upfront analysis.

---

## References

- [Stateless Agent Methodology](./stateless-agent-methodology.md) - This repository
- [Ralph Loop Orchestrator AGENTS.md](https://github.com/mikeyobrien/ralph-orchestrator/blob/main/AGENTS.md) - Source comparison document
