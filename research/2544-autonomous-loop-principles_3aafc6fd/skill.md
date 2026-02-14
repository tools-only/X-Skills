# Autonomous Loop Principles

Universal patterns for autonomous development systems, extracted from cross-framework analysis of 6 AI development frameworks during the ARL expert panel process (2026-02-12 to 2026-02-13).

---

## Overview

These 7 principles apply to any autonomous development loop — not just the ARL. They were discovered through comparative analysis of BMAD-METHOD, Gastown, GSD, Octocode, Ralph, and SAM. Each principle is supported by evidence from multiple frameworks with file:line citations.

---

## Scope

These principles complement the Stateless Agent Methodology (SAM) for iterative refinement scenarios. SAM handles single-pass pipelines (discover → plan → execute → verify). These principles address what happens when the pipeline must repeat: convergence detection, drift prevention, loop health monitoring.

**Reference:** [Stateless Agent Methodology](./stateless-software-engineering-framework.md) for SAM's canonical specification.

SAM is predominantly single-pass. When refinement loops are necessary, these principles provide the missing iteration-aware controls.

---

## Source Attribution

Extracted from cross-framework analysis using primary evidence from 6 repository source files with file:line citations. Method: multi-agent expert panel with cross-examination.

**Methodology:** [Expert Panel Methodology](./expert-panel-methodology.md)

**Evidence Source:** [ARL Expert Panel Q&A](../plugins/plugin-creator/skills/assessor/references/ARL/qa-expert-panel.md)

**Frameworks Analyzed:**

- BMAD-METHOD (8-step epic workflow)
- Gastown (molecule orchestration system)
- GSD (3-scale: project/phase/task)
- Octocode (GAN-inspired adversarial development)
- Ralph (event-driven hat routing)
- SAM (stateless agent methodology)

---

## Principle 1: Structure Over Instruction

**Statement:** Telling an AI agent to perform a check is unreliable. Structuring the pipeline so the check is the only possible path is reliable. Behavioral instructions cannot override architectural limitations.

**Evidence:**

All 6 frameworks converge on this principle through different structural mechanisms:

| Framework | Structural Enforcement Mechanism | Citation |
|-----------|----------------------------------|----------|
| SAM | Artifact-based pipeline stages | ssf:108 "Behavioral instructions cannot override architectural limitations" |
| BMAD | Sequential document chain | Each workflow step requires previous step's artifact |
| GSD | Context budget limits | Agents cannot exceed allocated context, preventing scope creep |
| Ralph | Backpressure gates | Machine-verifiable evidence (tests, lint) gates progression |
| Octocode | Adversarial flow | Generator/Discriminator tension forces quality structurally |
| Gastown | Zero Decisions in Code | Go code never decides, AI agents decide within structural constraints |

**Implication for autonomous loops:**

Gates in an autonomous loop must be structural (pipeline forces the check) rather than instructional (agent is told to check itself). Self-assessment without structural enforcement is unreliable.

**Counter-evidence:**

Octocode's STOP gates and verbalization requirements are behavioral instructions. This creates a tension: does Octocode work because of phase-gate structure (structural) or because of explicit STOP instructions (behavioral)? The expert panel did not resolve this tension.

---

## Principle 2: Front-Loading Reduces Runtime Human Gates

**Statement:** The more context captured before autonomous execution begins, the fewer human interventions required during execution. Front-loading converts runtime judgment calls into upfront decisions.

**Evidence:**

The 6 frameworks form a spectrum from deep to light front-loading, with inverse correlation to runtime human gates:

| Framework | Front-Loading Depth | Runtime Human Gates | Citation |
|-----------|---------------------|---------------------|----------|
| SAM | Deepest (4 stages before execution) | Single human touch at Discovery | ssf stages 1-4 |
| BMAD | Deep (sequential artifact chain) | 8 workflow steps with gates, each lighter due to captured context | workflow.xml |
| GSD | Medium (3-scale: project/phase/task) | 90% verify-type checkpoints (low judgment) | discuss-phase.md front-loaded decisions |
| Ralph | Medium (objective + specs) | 4 interaction points, objective-setting critical | hatless_ralph.rs:21-23 |
| Octocode | Light-medium (context.md + hints) | Fast-path can skip approval when 4 binary criteria pass | research:237-259 |
| Gastown | Structural (DAG at authoring time) | Runtime contact determined by role assignment | molecule.go:55-66 |

**What fails without front-loading:**

Each framework documents specific failure modes when front-loading is skipped:

- **SAM:** Agents hallucinate-fill missing information
- **GSD:** Technically valid but wrong-direction work
- **Octocode:** Review without understanding architecture
- **Gastown:** Escalate preventable issues

**Universal gap:**

No framework addresses re-triggering front-loading from within an iterative loop. All front-loading is one-shot. If execution reveals the original assessment was insufficient, no structured "return to front-loading" path exists.

---

## Principle 3: AI Cannot Reliably Self-Evaluate

**Statement:** An AI agent evaluating its own work shares the same blind spots that produced the work. Independent verification — structural separation between producer and evaluator — is necessary for reliable quality assessment.

**Evidence:**

| Framework | Independent Verification Mechanism | Citation |
|-----------|-----------------------------------|----------|
| SAM | Forensic Review Agent structurally separate from Execution Agent | ssf:856 "Self-confirmation bias → Independent verification agent" |
| Octocode | GAN-inspired adversarial flow: Generator vs Discriminator | MANIFEST.md:122-125 Cross-model validation eliminates shared biases |
| Ralph | Three-layer chain: Builder → Confessor → Handler | Handler spot-checks Confessor for fabricated issues |
| BMAD | Adversarial code review cross-references story claims vs git reality | code-review/instructions.xml:26-31 |
| Gastown | Witness zombie detection: agent self-report vs tmux session reality | Three agents cross-reference claims vs git state |

**Adversarial strength spectrum:**

From "check correctness" to "actively try to break":

| Strength Level | Framework | Mechanism | False Positive Risk |
|----------------|-----------|-----------|---------------------|
| Check correctness | SAM | Independent but not adversarial | Low |
| Find issues | Ralph | Confessor seeks problems | Medium |
| Forced minimum findings | BMAD | Must find 3-10 issues or halt | High |
| Zero-sum framing | Octocode | GAN-inspired: generator vs discriminator | Medium-High |

**Universal gap:**

No framework filters the reviewer's own false positives. All adversarial patterns are unidirectional. The producer cannot dispute a reviewer finding with counter-evidence.

---

## Principle 4: Compression Is Architectural, Not Instructional

**Statement:** The most effective information compression for human consumption occurs when the architecture forces compression — not when agents are instructed to summarize.

**Evidence:**

| Framework | Architectural Compression Mechanism | Citation |
|-----------|-------------------------------------|----------|
| GSD | Orchestrator operates at 10-15% context budget | execute-phase.md:322 physically cannot hold raw outputs |
| Gastown | Dashboard renders fixed Go struct (DashboardSummary) | HasAlerts boolean compresses all system state |
| Ralph | Check-in sends 4 scalar values: iteration, time, hat, cost | Telegram channel forces brevity |

**Contrast with instruction-based compression:**

| Framework | Instructional Compression (Less Effective) | Citation |
|-----------|-------------------------------------------|----------|
| Ralph | human.interact sends raw agent text | No compression on failure path |
| SAM | Artifacts self-contained with no cross-review aggregation | Human must read all artifacts |

**Implication for autonomous loops:**

A lead agent synthesizing multiple agents' outputs produces better compression when structurally constrained (limited context, fixed output format) rather than instructionally guided ("please summarize").

---

## Principle 5: Autonomous Loop Control Requires Iteration-Aware State

**Statement:** Single-pass frameworks cannot control iterative loops. The capabilities needed for autonomous loop operation (convergence tracking, oscillation detection, purpose drift detection) all require state that persists across iterations.

**Evidence from universal gaps:**

All 6 frameworks lack iteration-aware state because they are predominantly single-pass:

| Required Capability | Gap Across All Frameworks | Citation |
|---------------------|---------------------------|----------|
| Loop detection | Only Ralph has thrashing detection (count-based, 3 failures = escalate) | No framework detects output-similarity oscillation |
| Convergence tracking | No framework tracks finding-count-per-iteration | GSD tracks score delta within single phase only (qa-expert-panel Q3) |
| Purpose drift | No framework detects gradual semantic drift across iterations | Individually-acceptable changes cumulatively shift purpose |

**Why the gap exists:**

5 of 6 frameworks are predominantly single-pass pipelines (SAM, BMAD, Octocode) or bounded iterations with fresh context (GSD, Ralph). Requirements that emerge from iterative refinement (convergence, oscillation, cumulative drift) are invisible to single-pass designs.

**Implication:**

Any autonomous iterative refinement system must maintain cross-iteration state — a capability that must be designed from scratch because no surveyed framework provides it.

---

## Principle 6: Parallelism Enables Independent Verification Without Sequential Bottlenecks

**Statement:** Multiple independent agents checking different quality dimensions simultaneously provides stronger assurance than one agent checking dimensions sequentially, because independence reduces shared blind spots and parallelism reduces latency.

**Evidence:**

All 6 frameworks have sequential gates where serialization is an artifact of single-agent execution, not logical dependency. The expert panel identified 12+ parallelizable gates across frameworks.

| Framework | Parallel Pattern | Citation |
|-----------|------------------|----------|
| SAM | Spawn multiple focused workers, leader synthesizes | ssf:1134 explicit parallel specialist pattern |
| GSD | Spawn one debug agent per gap | diagnose-issues.md:87-88 |
| Octocode | Parallel domain writers with exclusive file ownership | doc-writer:50-56 |
| Gastown | Convoy formula spawns parallel polecats for different design dimensions | design.formula.toml:26-29 |

**Most common parallelizable gate type:**

"Multiple independent quality dimensions checked sequentially by one agent" — where dimensions operate on different aspects of the same artifact with no interdependency.

**Counter-argument (SAM):**

"Multi-agent microservices are an attractive trap because non-deterministic agents amplify coordination complexity; a monolithic loop can be more robust" (ssf:1070).

**Unresolved tension:**

Independence (reduces bias) vs coordination complexity (introduces failure modes). No framework has empirically measured the trade-off.

---

## Principle 7: Escalation Failure Paths Need More Compression Than Success Paths

**Statement:** All frameworks compress success paths well (structured reports, scores, tables). Failure paths — when the human most needs compressed context — are less compressed.

**Evidence:**

| Framework | Success Path Compression | Failure Path Compression | Citation |
|-----------|--------------------------|--------------------------|----------|
| Ralph | Structured check-in: 4 scalar values | Raw agent text on escalate.human | mod.rs:1958-2050 (no compression) |
| SAM | Structured artifacts with schemas | No severity classification — all NEEDS_WORK treated equally | Universal treatment |
| BMAD | Structured review reports | No mid-epic pattern detection — same finding repeated across stories | Retrospective aggregates late |
| GSD | Score + what's missing + next command | Best compressed failure path — diagnosis aggregation into table | Gap found presentation |

**Gastown's unique contribution:**

Time-based auto-escalation: Unacknowledged issues automatically become more urgent via stale_threshold (4h) with max_reescalations (2). Addresses the "forgotten escalation" failure mode — an issue that falls through the cracks is automatically re-surfaced with higher severity.

**No other framework implements this pattern.**

---

## Decision Tree: When Can a Human Gate Be Replaced?

A human gate can potentially be replaced by machine-verifiable conditions when ALL of the following hold:

1. **The check has an external truth source** — something outside the AI's own output to verify against (git state, test results, file existence, structural inventory).

   **Evidence:** BMAD git reality check (code-review/instructions.xml), Gastown tmux verification, Ralph backpressure evidence (event_parser.rs)

2. **The check operates on a single dimension** — one aspect of quality, not holistic judgment. Multi-dimensional judgment resists decomposition.

   **Evidence:** Parallelizable gates are single-dimension (qa-expert-panel Q5); irreducible gates are multi-dimensional

3. **Success and failure are distinguishable without domain knowledge** — the check can determine pass/fail from structural properties rather than semantic understanding.

   **Evidence:** Eliminable checks are structural (qa-expert-panel Q1); irreducible checks require domain judgment

4. **The scope is bounded** — the check operates on a defined artifact with known boundaries, not an open-ended assessment.

   **Evidence:** SAM's RT-ICA bounds scope before execution (ssf), GSD's context budget bounds plan scope (gsd-planner.md)

**When ANY of these conditions fails:**

The gate requires either human judgment or a more sophisticated verification mechanism (adversarial review, cross-examination between independent agents, or escalation).

---

## Unresolved Tensions

These tensions emerged from the expert panel discussion and remain unresolved:

### 1. Behavioral vs Structural Enforcement for Plan Approval

Octocode classifies plan approval as "reducible but not eliminable" — if plan quality were self-certifiable with evidence-backed traceability, approval could be auto-granted. No other framework claims plan approval is reducible.

**Source:** qa-expert-panel Q1 cross-examination note

**Resolution required:** Can machine-verified plan quality substitute for human plan approval?

### 2. Independence vs Coordination Complexity

More independent agents provide stronger verification but introduce coordination failure modes. SAM warns against "multi-agent microservices" while simultaneously documenting parallel specialist patterns.

**Source:** qa-expert-panel Q5 (SAM ssf:1070 vs ssf:1134)

**Resolution required:** What is the optimal balance? Under what conditions does independence outweigh coordination cost?

### 3. Adversarial Strength vs False Positive Rate

Stronger adversarial mandates (BMAD: "find 3-10 issues minimum") catch more real issues but also generate more false positives. No framework has found the right balance — none even measure the false positive rate.

**Source:** qa-expert-panel Q5 adversarial strength spectrum

**Resolution required:** What is the acceptable false positive rate? How to measure it?

### 4. Front-Loading Completeness vs Front-Loading Cost

SAM front-loads the most (4 stages) but this delays execution. Ralph front-loads the least (objective only) but requires more runtime intervention. The optimal amount of front-loading depends on task characteristics that no framework formally classifies.

**Source:** qa-expert-panel Q2 front-loading spectrum

**Resolution required:** What task characteristics determine optimal front-loading depth? How to classify tasks for this purpose?

---

## References

**Primary Evidence Source:**

- [ARL Expert Panel Q&A](../plugins/plugin-creator/skills/assessor/references/ARL/qa-expert-panel.md) — Cross-framework analysis with file:line citations

**Complementary Analysis:**

- [General Theory Synthesis](../plugins/plugin-creator/skills/assessor/references/ARL/synthesis-general-theory.md) — Extended principle descriptions

**Methodology:**

- [Expert Panel Methodology](./expert-panel-methodology.md) — Multi-agent cross-examination process

**Framework Integration:**

- [Stateless Agent Methodology](./stateless-software-engineering-framework.md) — SAM's canonical specification
