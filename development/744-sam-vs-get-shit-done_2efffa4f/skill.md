# Stateless Agent Methodology vs Get Shit Done (GSD)

A comparative analysis of two approaches to working effectively with LLM agents.

---

## Executive Summary

Both methodologies address the same core problem: **LLM agents optimize for apparent completion over correct completion**, leading to quality degradation during extended work sessions.

| Aspect       | Stateless Agent Methodology                       | Get Shit Done (GSD)                                   |
| ------------ | ------------------------------------------------- | ----------------------------------------------------- |
| **Origin**   | Theoretical framework from observed failure modes | Production system (npm: `get-shit-done-cc`)           |
| **Target**   | Any LLM agent workflow                            | Claude Code specifically                              |
| **Maturity** | Conceptual methodology                            | v1.9.13 (local checkout) with 27 commands + 11 agents |
| **Focus**    | Why agents fail and how to structure work         | Executable workflow with automation                   |

---

## Problem Statement Alignment

Both methodologies identify the same fundamental issues:

### Shared Problem Analysis

| Problem                                    | Stateless Agent Methodology                                                                                                               | GSD                                                                                                     |
| ------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------- |
| Long-context degradation (“context rot”)   | Performance can degrade as context length increases (including “lost in the middle” effects); mitigate by bounding context per task/stage | “Context rot” — quality degrades as context fills; mitigate via small plans and fresh executor contexts |
| Training data staleness (knowledge cutoff) | One-shot approaches drift into stale priors unless grounded/verified                                                                      | Research phase investigates domain before planning; fresh executor context reduces accumulated drift    |
| Knowledge gaps                             | Work involves recent or internal knowledge not in training                                                                                | Research phase investigates domain before planning                                                      |
| Apparent vs actual completion              | "Disables tests to achieve 'all tests passing'"                                                                                           | Verification phase independently confirms goals                                                         |

### The Core Dysfunction

**Stateless Agent Methodology**:

> "Claude is not a knowledge worker—Claude is a stateless computation engine. LLM agents cannot perform JIT identification of knowledge gaps through self-reflection."

**GSD**:

> "Context engineering hides complexity — plans capped at 2-3 tasks, fresh context per plan."

Both recognize that agents cannot self-assess, so the solution must externalize assessment.

---

## Architectural Comparison

### Phase Structure

| Stateless Agent Methodology  | GSD Equivalent                                          |
| ---------------------------- | ------------------------------------------------------- |
| Stage 1. Discovery           | `/gsd:discuss-phase` → creates `{phase}-CONTEXT.md`     |
| Stage 2. Planning (RT-ICA)   | `/gsd:plan-phase` with research agents                  |
| Stage 3. Context Integration | Codebase mapper agent + architecture research           |
| Stage 4. Task Decomposition  | Planner agent creates XML plan structure                |
| Stage 5. Execution           | `/gsd:execute-phase` with wave-based parallel execution |
| Stage 6. Forensic Review     | `/gsd:verify-work` + verifier agent                     |
| Stage 7. Orchestration Loop  | Built into workflow orchestration                       |
| Stage 8. Final Verification  | `/gsd:complete-milestone` with goal verification        |

### Key Structural Differences

| Dimension             | Stateless Agent Methodology               | GSD                                              |
| --------------------- | ----------------------------------------- | ------------------------------------------------ |
| **Abstraction level** | Methodology (principles)                  | Implementation (commands + agents)               |
| **Execution model**   | Sequential stages, fresh sessions         | Wave-based parallel execution                    |
| **State management**  | Implicit (task files contain all context) | Explicit (`STATE.md` tracks position, decisions) |
| **Granularity**       | 7 stages + orchestration loop             | 6-step cycle × N phases per milestone            |
| **Human interaction** | Stage boundaries                          | Checkpoints within stages                        |

---

## Detailed Phase Mapping

### Phase 1: Discovery → Discuss Phase

**Stateless Agent Methodology**:

- User discussion following data-gathering framework
- Output: Feature requirements guide, NFR guide, notes

**GSD**:

- `/gsd:discuss-phase N` captures implementation preferences
- Identifies gray areas (visual features, API format, error handling)
- Output: `{phase}-CONTEXT.md`

**Comparison**: GSD is more structured — focuses specifically on gray areas and preferences that would otherwise become generic defaults. Stateless methodology is broader but less prescriptive.

---

### Phase 2: Assessment → Research + Planning

**Stateless Agent Methodology**:

- Review goals, identify prerequisites
- Block planning until prerequisites verified
- Output: Feature design guide with verified prerequisites

**GSD**:

- 4 parallel research agents: stack, feature, architecture, pitfalls
- Creates `{phase}-RESEARCH.md`
- Plan verification before execution
- Output: Multiple `{phase}-{N}-PLAN.md` files

**Comparison**: GSD separates research (investigation) from planning (task creation). Stateless methodology combines prerequisite verification into a single assessment phase. GSD's parallel research is more sophisticated; stateless methodology's blocking requirement is more explicit.

---

### Phase 3: Integration → Codebase Mapping

**Stateless Agent Methodology**:

- Map plan against existing code
- Identify conflicts, find reusable systems
- Add file/URL references
- Output: Contextualized implementation plan

**GSD**:

- `gsd-codebase-mapper` agent analyzes existing code
- Architecture research identifies patterns
- Plan files include `<files>` section with paths
- Output: Plans grounded in actual codebase

**Comparison**: Both require grounding plans in reality. GSD integrates this into the research phase; stateless methodology treats it as a distinct step. The intent is identical.

---

### Phase 4: Task Generation → Plan Creation

**Stateless Agent Methodology**:

- Decompose into discrete tasks
- Follow: interfaces → failing tests → (implement → lint → test)
- Output: Task files with explicit constraints

**GSD**:

- Planner agent creates 2-5 atomic plans per phase
- XML structure: `<name>`, `<files>`, `<action>`, `<verify>`, `<done>`
- Each plan: 2-3 tasks maximum
- Output: `{phase}-{N}-PLAN.md` files

**Comparison**: GSD enforces small plan size (2-3 tasks) to prevent context degradation. Stateless methodology describes the pattern but doesn't enforce size constraints. GSD's XML format is concrete; stateless methodology is format-agnostic.

---

### Phase 5: Execution → Execute Phase

**Stateless Agent Methodology**:

- Agent receives task file as complete prompt
- Task contains all necessary information
- Self-verification steps throughout
- Output: Completed work + self-verification

**GSD**:

- Fresh, bounded context per executor
- Wave-based parallel execution (pre-computed dependencies)
- One atomic commit per task
- Output: `{phase}-{N}-SUMMARY.md` + commits

**Comparison**: GSD's wave execution enables parallelism; stateless methodology is sequential. Both emphasize complete context per execution unit. GSD adds atomic commits as a key output; stateless methodology doesn't address version control.

---

### Phase 6: Forensics → Verify Work

**Stateless Agent Methodology**:

- Independent verification of task completion
- Validate against plan
- Report status objectively
- Output: Complete/incomplete determination

**GSD**:

- `/gsd:verify-work N` runs UAT
- Extract testable deliverables
- Auto-diagnose failures with debugger agents
- Output: `{phase}-UAT.md` + fix plans if needed

**Comparison**: GSD includes automated debugging and fix generation. Stateless methodology stops at status determination. GSD's gap closure loop (`--gaps` flag) handles incomplete work automatically.

---

### Phase 7-8: Orchestration + Verification → Complete Milestone

**Stateless Agent Methodology**:

- Orchestrator creates follow-up tasks
- Dispatches workers to incomplete tasks
- Final verification against goal and acceptance criteria
- Output: Feature sign-off or additional tasks

**GSD**:

- Workflow orchestration coordinates agents throughout
- `/gsd:complete-milestone` archives and tags
- `/gsd:new-milestone` starts fresh cycle
- Output: Git tag, archived milestone, next milestone setup

**Comparison**: GSD treats orchestration as continuous (throughout execution), not a separate phase. Milestone concept provides natural boundaries for releases.

---

## Key Differentiators

### What GSD Has That Stateless Methodology Lacks

| Feature                 | GSD Implementation                     | Stateless Gap                               |
| ----------------------- | -------------------------------------- | ------------------------------------------- |
| **Executable commands** | 27 slash commands (local checkout)     | No implementation                           |
| **Specialized agents**  | 11 agents with specific roles          | Agent separation described, not implemented |
| **Wave parallelism**    | Pre-computed dependency execution      | Sequential only                             |
| **Git integration**     | Atomic commits per task                | Not addressed                               |
| **State persistence**   | `STATE.md` tracks across sessions      | Relies on task files only                   |
| **Gap closure loop**    | Automatic plan generation for failures | Manual intervention required                |
| **Configuration**       | Model profiles, workflow settings      | No configuration model                      |
| **Templates**           | 20+ structural templates               | No output formats defined                   |

### What Stateless Methodology Emphasizes That GSD Assumes

| Concept                       | Stateless Emphasis             | GSD Approach                               |
| ----------------------------- | ------------------------------ | ------------------------------------------ |
| **Why agents fail**           | Detailed failure mode analysis | Assumes understanding, focuses on solution |
| **Self-assessment inability** | Core thesis, explicit          | Implicit in design                         |
| **Training data distrust**    | Central problem statement      | Addressed via fresh context, not discussed |
| **Methodology vs process**    | Principles first               | Process first                              |
| **Agent-agnostic**            | Any LLM agent                  | Claude Code specific                       |

---

## Failure Mode Coverage

| Failure Mode                               | Stateless Mitigation                                          | GSD Mitigation                                                 |
| ------------------------------------------ | ------------------------------------------------------------- | -------------------------------------------------------------- |
| Training data hallucination / stale priors | Grounding + explicit constraints + deterministic backpressure | Research phase + deterministic checks + fresh executor context |
| Skipping prerequisites                     | Assessment phase blocks planning                              | Discuss phase captures preferences                             |
| Apparent vs actual completion              | Forensic phase verification                                   | Verifier agent + UAT                                           |
| Long-context degradation (“context rot”)   | Bounded context per task/stage + artifacts as source of truth | Fresh context per plan (smaller unit) + small plan sizes       |
| Rationalizing out of process               | Process is the task                                           | Process is automated (commands)                                |

---

## Complementary Strengths

### Stateless Methodology Contributes

1. **Theoretical foundation** — Explains _why_ the approach works
2. **Agent-agnostic principles** — Applicable beyond Claude Code
3. **Explicit failure mode analysis** — Diagnostic framework for new problems
4. **Self-assessment thesis** — Core insight that drives design decisions

### GSD Contributes

1. **Production-ready implementation** — Immediately usable
2. **Parallel execution model** — Performance optimization
3. **Git integration** — Version control as first-class concern
4. **Debugging automation** — Self-healing workflow
5. **Configuration system** — Adaptable to different needs

---

## Synthesis: A Unified View

Both systems share the same architecture when viewed abstractly:

```
┌─────────────────────────────────────────────────────────────────┐
│                    STATELESS AGENT PATTERN                      │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  ┌──────────┐    ┌──────────┐    ┌──────────┐    ┌──────────┐  │
│  │ DISCOVER │───▶│  PLAN    │───▶│ EXECUTE  │───▶│ VERIFY   │  │
│  │          │    │          │    │          │    │          │  │
│  │ Complete │    │ Complete │    │ Complete │    │ Complete │  │
│  │ context  │    │ context  │    │ context  │    │ context  │  │
│  └──────────┘    └──────────┘    └──────────┘    └──────────┘  │
│       │               │               │               │         │
│       ▼               ▼               ▼               ▼         │
│  ┌──────────┐    ┌──────────┐    ┌──────────┐    ┌──────────┐  │
│  │ Artifact │    │ Artifact │    │ Artifact │    │ Artifact │  │
│  │ (reqs)   │    │ (plans)  │    │ (code)   │    │ (report) │  │
│  └──────────┘    └──────────┘    └──────────┘    └──────────┘  │
│                                                                 │
│  Key insight: Each phase receives COMPLETE context.             │
│  No phase depends on agent "remembering" anything.              │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

**The pattern**:

1. Each phase gets complete context (stateless)
2. Each phase produces persistent artifacts
3. Artifacts become input to next phase
4. Assessment is externalized (separate agents/phases)
5. Verification is independent from execution

---

## Recommendations

### If Using Stateless Agent Methodology

Consider adopting from GSD:

- Wave-based parallelism for execution
- Atomic commits per task
- Gap closure loop for incomplete work
- Explicit state file across sessions
- Size constraints on task units (2-3 tasks max)

### If Using GSD

The stateless methodology provides:

- Theoretical justification for why GSD works
- Framework for extending to non-Claude agents
- Failure mode vocabulary for debugging new issues
- Principles for designing new workflows

### For New System Design

Combine:

1. **Stateless methodology's principles** — Design rationale
2. **GSD's structure** — Milestone/phase/plan hierarchy
3. **GSD's execution model** — Wave parallelism
4. **Stateless methodology's forensics** — Independent verification
5. **GSD's automation** — Debugging and gap closure

---

## Conclusion

The Stateless Agent Methodology and Get Shit Done are **the same insight at different abstraction levels**.

- **Stateless Agent Methodology** = Why this works (theory)
- **Get Shit Done** = How to do it (practice)

Neither is complete without the other:

- Theory without implementation is academic
- Implementation without theory is fragile

Together they form a comprehensive approach to reliable LLM agent workflows.

---

## References

- [get-shit-done-cc (npm)](https://www.npmjs.com/package/get-shit-done-cc) (accessed 2026-01-26)
- [get-shit-done-cc versions (npm)](https://www.npmjs.com/package/get-shit-done-cc?activeTab=versions) (accessed 2026-01-26)
- [glittercowboy/get-shit-done (GitHub)](https://github.com/glittercowboy/get-shit-done) (accessed 2026-01-26)
- Commands list (GitHub): [commands/gsd/](https://github.com/glittercowboy/get-shit-done/tree/main/commands/gsd) (accessed 2026-01-26)
- Agents list (GitHub): [agents/](https://github.com/glittercowboy/get-shit-done/tree/main/agents) (accessed 2026-01-26)
- Repo evidence (accessed 2026-01-26): [`package.json`](https://github.com/glittercowboy/get-shit-done/blob/main/package.json) (version `1.9.13`), plus counts observed in `commands/gsd/` (27) and `agents/` (11) from a local checkout of the same repository.
