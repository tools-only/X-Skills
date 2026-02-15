---
description: Autonomous Refinement Loop for iterative skill improvement — Assess, Plan, Implement, Review, Repeat until converged. Machine-verifiable gates (R1-R10) replace human judgment at each phase. Use when refining skills through multiple autonomous iterations.
user-invocable: true
model: opus
context: fork
---

# Autonomous Refinement Loop (ARL)

The ARL is a logical process for iterative skill refinement that operates without human intervention by replacing judgment gates with machine-verifiable conditions.

**Core loop:** Assess → Plan → Implement → Review → Repeat until converged or stopped.

## When to Use

Activate the ARL when:

- A skill has quality issues requiring multiple improvement passes
- Iterative refinement is needed with measurable progress tracking
- Human-out-of-loop operation is desired for bounded improvement tasks
- Assessor has identified findings that require multiple iterations to resolve
- Convergence must be tracked objectively across iterations

**Prerequisites before activation:** Sufficient information must be available upfront to operate the loop without runtime escalation for missing context (R1 verification required).

## Prerequisites: Information Completeness Gate (R1)

Before entering the loop, verify information completeness using the RT-ICA pattern (Reverse Thinking - Information Completeness Assessment). The loop cannot start without sufficient input.

**What must be present:**

- Skill's stated purpose and intended use cases
- Initial assessment findings with file:line evidence
- Refinement plan with acceptance criteria
- Knowledge of downstream references (what depends on this skill)
- Quality baselines (structural inventory before changes)

**What R1 checks:** Classify each information need as available, derivable, or missing. Block progression if critical information is missing.

**When R1 activates:**

- At loop entry (iteration 0)
- At re-entry after escalation

**Failure state:** Loop begins with missing information, produces wrong output that must be reverted, or escalates for information that should have been captured upfront.

**SOURCE:** [Synthesis: ARL-Applicable](./references/synthesis-arl-applicable.md) lines 11-32, [Human-Out-of-Loop Prerequisites](./references/human-out-of-loop-prerequisites.md) lines 166-177

## Loop Structure

One ARL iteration follows this gate activation sequence:

```
R1 (info completeness) → R5 (purpose anchor) → ASSESS → R3 (validity filtering) →
R7+R2 (convergence/loop check) → PLAN → R4+R8 (quality/proportionality) →
R10 (split justification) → IMPLEMENT → R6 (content-loss) → R9 (downstream impact) → Loop
```

**Gate timing:**

- **Before assessment:** R1 (information sufficient?), R5 (purpose drift?)
- **After assessment:** R3 (findings valid?), R7 (converging?), R2 (oscillating/stalled?)
- **After planning:** R4 (plan sound?), R8 (changes proportional?), R10 (splits justified?)
- **After implementation:** R6 (content preserved?), R9 (references still valid?)

**SOURCE:** [Synthesis: ARL-Applicable](./references/synthesis-arl-applicable.md) lines 260-305

## The 10 Gates

| Gate | What It Checks | When It Fires | What Failure Looks Like |
|------|----------------|---------------|------------------------|
| **R1: Information Completeness** | Sufficient context to operate loop without escalation | Loop entry, re-entry after escalation | Loop proceeds with gaps, agent hallucinate-fills missing information, produces fluent but wrong artifacts |
| **R2: Loop Detection** | Oscillating, stalling, or exceeding resource bounds | Start of each iteration before assessment | Loop runs indefinitely without converging, fix A breaks B repeatedly, same findings recurring |
| **R3: Validity Filtering** | Findings have verifiable evidence (file:line citations) | After assessment, before planning | False positives consume iteration budget, regressions introduced, phantom issues trigger changes |
| **R4: Plan Quality** | Plan internally consistent, addresses actual findings | After planning, before implementation | Inconsistent plan proceeds, addresses wrong findings, changes must be reverted |
| **R5: Purpose Anchor** | Skill still serves original stated purpose | Captured at iteration 0, checked each iteration | After N iterations, skill optimized for assessor metrics but no longer serves original use case |
| **R6: Content-Loss Detection** | All semantic units preserved after changes | After implementation, before next iteration | Refactoring removes sections deemed "redundant", no gate catches removal, human discovers loss later |
| **R7: Convergence Tracking** | Findings decreasing, stable, or alternating across iterations | Each iteration boundary after assessment | Loop cannot determine progress, fixes trivial issues indefinitely, or oscillates without converging |
| **R8: Proportionality Check** | Proposed fix proportional to finding severity | During plan quality gate (R4) | Low-severity finding triggers high-scope change that introduces risk without proportional benefit |
| **R9: Downstream Impact** | All references still resolve after changes | After implementation, alongside R6 | Refactoring renames file, breaks three other skills linking to old path, not detected until invoked |
| **R10: Split Justification** | New skill independently viable, not just parent-dependent | When plan proposes splitting content | Skill split into three pieces, two only invoked from parent, adds navigation complexity without value |

**SOURCE:** [Synthesis: ARL-Applicable](./references/synthesis-arl-applicable.md) — each gate documented in sections R1-R10

## Exit Conditions

The loop stops when:

- **Convergence (R7):** Finding count reaches zero or stabilizes below value threshold
- **Failure state (R2):** Oscillation detected, no progress across N iterations
- **Insufficient info (R1):** Cannot proceed without human input
- **Max iterations:** Safety limit reached (prevent unbounded execution)

**Escalation vs termination:** R2 failure states may trigger self-correction if recoverable, or escalate with diagnosis if not. R1 gaps always escalate to request missing context.

**SOURCE:** [Synthesis: ARL-Applicable](./references/synthesis-arl-applicable.md) lines 304-305

## Framework Patterns

The 10 R-requirements map to existing framework mechanisms at three coverage levels:

| Coverage | R-Requirements | Framework Patterns Available |
|----------|---------------|------------------------------|
| **Import directly** | R1, R3, R4 | RT-ICA (SAM), GAN-inspired validation (Octocode/BMAD), 7-dimension plan checking (GSD) |
| **Partial coverage** | R2, R5, R9 | Bounded iteration count (GSD), objective injection (Ralph), downstream impact analysis (Octocode) |
| **Build from scratch** | R6, R7, R8, R10 | No framework provides content-loss detection, convergence tracking, proportionality checks, or split justification |

**Key insights:**

- **R1 (RT-ICA):** SAM's RT-ICA is one-shot; ARL needs re-triggerable RT-ICA within the loop
- **R2 (Loop detection):** No framework detects finding-level oscillation across iterations
- **R5 (Purpose anchor):** GSD has deviation rules within execution; ARL needs cumulative drift detection across iterations
- **R6 (Content-loss):** Entirely novel — no framework compares structural inventories before/after changes
- **R7 (Convergence):** Entirely novel — requires cross-iteration state tracking absent from single-pass pipelines
- **R8 (Proportionality):** Novel comparison of finding severity vs change scope
- **R10 (Split justification):** Novel independent viability assessment for extracted skills

**SOURCE:** [Synthesis: ARL-Applicable](./references/synthesis-arl-applicable.md) — each R-requirement section includes framework pattern table

## Universal Principles

The ARL applies seven universal patterns for autonomous development systems:

1. **Structure Over Instruction:** Pipeline forces checks rather than asking agents to check themselves
2. **Front-Loading Reduces Runtime Gates:** More context captured upfront = fewer human interventions during execution
3. **AI Cannot Reliably Self-Evaluate:** Independent verification structurally separates producer from evaluator
4. **Compression Is Architectural:** Information compression works when architecture forces it, not when instructed
5. **Iteration-Aware State Required:** Loop control needs state persisting across iterations (convergence, oscillation, drift)
6. **Parallelism Enables Independent Verification:** Multiple agents checking different dimensions simultaneously reduces shared blind spots
7. **Failure Paths Need More Compression:** Escalation context needs better compression than success paths

**Cross-reference:** [Synthesis: General Theory](./references/synthesis-general-theory.md) documents these patterns with evidence from 6 surveyed frameworks

## References

Complete detail on each gate, framework patterns, and prerequisites:

- [Synthesis: ARL-Applicable](./references/synthesis-arl-applicable.md) — R1-R10 mapped to framework mechanisms, loop structure
- [Synthesis: General Theory](./references/synthesis-general-theory.md) — 7 universal principles for autonomous development systems
- [ARL Research](./references/autonomous-refinement-loop-research.md) — Problem definition, what replaces human judgment
- [Human-Out-of-Loop Prerequisites](./references/human-out-of-loop-prerequisites.md) — Conditions for removing human gates
- [Expert Panel Q&A](./references/qa-expert-panel.md) — Phase 1-2 expert panel findings
- [Autonomous Loop Principles](https://github.com/bitflight-devops/stateless-agent-methodology/blob/main/autonomous-loop-principles.md) — SAM extension analysis (cross-reference)
