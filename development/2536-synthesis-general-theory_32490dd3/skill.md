# Synthesis: General Theory — Universal Patterns for Autonomous Development Systems

This document synthesizes universal patterns discovered by the ARL expert panel that apply to any autonomous development loop — not just the ARL. Every claim traces to Phase 1 Q&A evidence and Phase 2 R-requirement mapping.

**Scope boundary:** This document describes what and why. It does not describe how to build it. No schemas, pseudocode, thresholds, or file paths.

---

## Principle 1: Structure Over Instruction

**Statement:** Telling an AI agent "please do X" is unreliable. Structuring the pipeline so X is the only possible path is reliable.

**Evidence:** All 6 surveyed frameworks converge on this principle through different mechanisms:

- SAM states it explicitly: "Behavioral instructions cannot override architectural limitations" (ssf:108)
- BMAD enforces via sequential document chain — each step requires the previous step's artifact
- GSD enforces via context budget — agents cannot exceed allocated context, preventing scope creep structurally
- Ralph enforces via backpressure gates — machine-verifiable evidence (tests, lint) gates progression
- Octocode enforces via adversarial flow — Generator/Discriminator tension forces quality structurally
- Gastown enforces via ZFC (Zero Decisions in Code) — Go code never decides, AI agents decide within structural constraints

**Source:** Q1 (convergence pattern across all 6 frameworks), Q3 (branching point analysis showing AI judgment dominance).

**Implication for autonomous loops:** Gates in an autonomous loop must be structural (pipeline forces the check) rather than instructional (agent is told to check itself). Self-assessment without structural enforcement is unreliable.

---

## Principle 2: Front-Loading Reduces Runtime Human Gates

**Statement:** The more context captured before autonomous execution begins, the fewer human interventions required during execution. Front-loading converts runtime judgment calls into upfront decisions.

**Evidence:** The 6 frameworks form a spectrum from deep to light front-loading, with inverse correlation to runtime human gates:

- SAM (deepest — 4 stages before execution): Human touches system once at Discovery, then pipeline runs autonomously
- BMAD (deep — sequential artifact chain): 8 workflow steps with human gates, but each gate is lighter because previous stages captured context
- GSD (3-scale — project/phase/task): 90% of runtime checkpoints are verify-type (low judgment) because discuss-phase front-loaded the decisions
- Ralph (medium — objective + specs): 4 interaction points, but objective-setting is the critical one; runtime gates handle exceptions
- Octocode (light-medium — context.md + hints): Fast-path can skip human approval entirely when 4 binary criteria pass
- Gastown (structural — DAG at authoring time): Runtime human contact determined by role assignment, not situation assessment

**Source:** Q2 (front-loading spectrum table), Q1 (interaction point inventory).

**What fails without front-loading:** Each framework documents specific failure modes when front-loading is skipped (Q2-Q4): agents hallucinate-fill missing information (SAM), produce technically valid but wrong-direction work (GSD), review without understanding architecture (Octocode), escalate preventable issues (Gastown).

**Gap in current theory:** No framework addresses re-triggering front-loading from within an iterative loop. All front-loading is one-shot. If execution reveals the original assessment was insufficient, no structured "return to front-loading" path exists.

---

## Principle 3: AI Cannot Reliably Self-Evaluate

**Statement:** An AI agent evaluating its own work shares the same blind spots that produced the work. Independent verification — structural separation between producer and evaluator — is necessary for reliable quality assessment.

**Evidence:**

- SAM: "Self-confirmation bias → Independent verification agent; do not rely on self-critique alone" (ssf:856). Forensic Review Agent is structurally separate from Execution Agent.
- Octocode: GAN-inspired adversarial flow explicitly separates Generator (produces) from Discriminator (evaluates). Cross-model validation eliminates shared model biases (MANIFEST.md:122-125).
- Ralph: Three-layer chain (Builder → Confessor → Handler) where each layer checks the previous. Handler spot-checks Confessor — if Confessor fabricated issues, Handler escalates.
- BMAD: Adversarial code review uses git reality as an external truth source — cross-references story claims against actual file state (code-review/instructions.xml:26-31).
- Gastown: Witness zombie detection cross-references agent self-reported state against tmux session reality. Three agents' claims cross-referenced against git state for merge verification.

**Source:** Q1 (root cause analysis — all frameworks identify AI cannot self-evaluate), Q5 (adversarial pattern comparison).

**Adversarial strength spectrum:** From "check correctness" (SAM — independent but not adversarial) to "actively try to break" (Octocode — GAN zero-sum framing, BMAD — forced minimum findings with halt on zero). Stronger adversarial mandates catch more real issues but also produce more false positives.

**Universal gap:** No framework filters the reviewer's own false positives. All adversarial patterns are unidirectional. The producer cannot dispute a reviewer finding with counter-evidence.

---

## Principle 4: Compression Is Architectural, Not Instructional

**Statement:** The most effective information compression for human consumption occurs when the architecture forces compression — not when agents are instructed to summarize.

**Evidence:**

- GSD's orchestrator operates at 10-15% context (execute-phase.md:322), which physically cannot hold raw agent outputs. Compression is not optional — it's a consequence of the architecture.
- Gastown's dashboard renders a fixed Go struct (DashboardSummary with HasAlerts boolean), compressing all system state into a predetermined format.
- Ralph's check-in sends 4 scalar values (iteration, time, hat, cost) — extreme compression because the Telegram channel forces brevity.

**Contrast with instruction-based compression:**

- Ralph's human.interact sends raw agent text with no compression — the agent is not structurally constrained.
- SAM's artifacts are self-contained with no cross-review aggregation — each artifact is complete but the human must read them all.

**Source:** Q5 (escalation compression comparison), Q2 (information presentation patterns).

**Implication for autonomous loops:** A lead agent synthesizing multiple agents' outputs produces better compression when structurally constrained (limited context, fixed output format) rather than instructionally guided ("please summarize").

---

## Principle 5: Autonomous Loop Control Requires Iteration-Aware State

**Statement:** Single-pass frameworks cannot control iterative loops. The capabilities needed for autonomous loop operation (convergence tracking, oscillation detection, purpose drift detection) all require state that persists across iterations.

**Evidence:**

- R2 (Loop detection): Only Ralph has any thrashing detection, and it's count-based (3 failures = escalate). No framework detects output-similarity oscillation across iterations.
- R7 (Convergence tracking): No framework tracks finding-count-per-iteration. GSD tracks score delta but within a single phase, not across loop iterations.
- R5 (Purpose drift): No framework detects gradual semantic drift where individually-acceptable changes cumulatively shift purpose across iterations.

**Source:** Q3 (universal gap analysis — R2, R7, R8, R10 absent from all frameworks), Q4 (drift detection universally absent).

**Why the gap exists:** 5 of 6 frameworks are predominantly single-pass pipelines (SAM, BMAD, Octocode) or bounded iterations with fresh context (GSD, Ralph). The requirements that emerge from iterative refinement (convergence, oscillation, cumulative drift) are invisible to single-pass designs because they never encounter them.

**Implication:** Any autonomous iterative refinement system must maintain cross-iteration state — a capability that must be designed from scratch because no surveyed framework provides it.

---

## Principle 6: Parallelism Enables Independent Verification Without Sequential Bottlenecks

**Statement:** Multiple independent agents checking different quality dimensions simultaneously provides stronger assurance than one agent checking dimensions sequentially, because independence reduces shared blind spots and parallelism reduces latency.

**Evidence:**

- All 6 frameworks have sequential gates where serialization is an artifact of single-agent execution, not logical dependency (Q5 Q2 analysis — 12+ parallelizable gates identified across frameworks).
- The most common parallelizable gate type: "multiple independent quality dimensions checked sequentially by one agent" — where dimensions operate on different aspects of the same artifact with no interdependency.
- SAM documents parallel specialist pattern explicitly: "spawn multiple focused workers, each produces a bounded artifact, then a leader synthesizes" (ssf:1134).
- GSD implements parallel diagnostic agents: "spawn one debug agent per gap" (diagnose-issues.md:87-88).
- Octocode implements parallel domain writers with exclusive file ownership (doc-writer:50-56).
- Gastown's convoy formula spawns parallel polecats for different design dimensions (design.formula.toml:26-29).

**Source:** Q5 (multi-agent pattern inventory, parallelization opportunities).

**SAM's counter-argument:** "Multi-agent 'microservices' are an attractive trap because non-deterministic agents amplify coordination complexity; a monolithic loop can be more robust" (ssf:1070). This tension between independence (reduces bias) and coordination complexity (introduces failure modes) is unresolved across the surveyed frameworks.

---

## Principle 7: Escalation Failure Paths Need More Compression Than Success Paths

**Statement:** All frameworks compress success paths well (structured reports, scores, tables). Failure paths — when the human most needs compressed context — are less compressed.

**Evidence:**

- Ralph sends raw agent text on escalate.human — no compression on the failure path (mod.rs:1958-2050).
- SAM has no severity classification — all NEEDS_WORK issues treated equally.
- BMAD has no mid-epic pattern detection — the human sees the same review finding repeated across stories before the retrospective aggregates it.
- GSD's gap_found presentation is the most compressed failure path — score + what's missing + next command. But diagnosis aggregation (parallel debug agents → table) is the strongest example.

**Source:** Q5 (escalation compression comparison table).

**Gastown's unique contribution — time-based auto-escalation:** Unacknowledged issues automatically become more urgent via stale_threshold (4h) with max_reescalations (2). This addresses the "forgotten escalation" failure mode — an issue that falls through the cracks is automatically re-surfaced with higher severity. No other framework implements this.

---

## Decision Tree: When Can a Human Gate Be Replaced?

Based on the cross-framework evidence, a human gate can potentially be replaced by machine-verifiable conditions when ALL of the following hold:

1. **The check has an external truth source** — something outside the AI's own output to verify against (git state, test results, file existence, structural inventory). [Q1: BMAD git reality check, Gastown tmux verification, Ralph backpressure evidence]

2. **The check operates on a single dimension** — one aspect of quality, not holistic judgment. Multi-dimensional judgment resists decomposition. [Q5: parallelizable gates are single-dimension; irreducible gates are multi-dimensional]

3. **Success and failure are distinguishable without domain knowledge** — the check can determine pass/fail from structural properties rather than semantic understanding. [Q1: eliminable checks are structural; irreducible checks require domain judgment]

4. **The scope is bounded** — the check operates on a defined artifact with known boundaries, not an open-ended assessment. [Q2: SAM's RT-ICA bounds scope before execution; GSD's context budget bounds plan scope]

When ANY of these conditions fails, the gate requires either human judgment or a more sophisticated verification mechanism (adversarial review, cross-examination between independent agents, or escalation).

**Source:** Q1 (eliminable/front-loadable/irreducible classification), Q3 (detection mechanism taxonomy), Q4 (machine-checkable vs judgment-dependent comparison).

---

## Unresolved Tensions

These tensions emerged from the expert panel discussion and remain unresolved:

1. **Behavioral vs structural enforcement for plan approval.** Octocode classifies plan approval as "reducible but not eliminable" — if plan quality were self-certifiable with evidence-backed traceability, approval could be auto-granted. No other framework claims plan approval is reducible. [Q1: cross-examination note]

2. **Independence vs coordination complexity.** More independent agents provide stronger verification but introduce coordination failure modes. SAM warns against "multi-agent microservices" while simultaneously documenting parallel specialist patterns. [Q5: SAM ssf:1070 vs ssf:1134]

3. **Adversarial strength vs false positive rate.** Stronger adversarial mandates (BMAD: "find 3-10 issues minimum") catch more real issues but also generate more false positives. No framework has found the right balance — none even measures the false positive rate. [Q5: adversarial strength spectrum]

4. **Front-loading completeness vs front-loading cost.** SAM front-loads the most (4 stages) but this delays execution. Ralph front-loads the least (objective only) but requires more runtime intervention. The optimal amount of front-loading depends on task characteristics that no framework formally classifies for this purpose. [Q2: front-loading spectrum]
