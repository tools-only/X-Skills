# Synthesis: ARL-Applicable — R1–R10 Mapped to Framework Mechanisms

This document maps each ARL requirement (R1–R10) to what the logical process should do, when it activates, why it's needed, what success/failure looks like, and which framework patterns inform it.

**Scope boundary:** This document describes what and why. It does not describe how to build it. No schemas, pseudocode, thresholds, file paths, or algorithm specifications.

**Source:** All claims trace to Phase 1 Q&A exchanges and Phase 2 R-requirement mapping in the Q&A file (qa-expert-panel.md).

---

## R1 — Information Completeness Gate (RT-ICA Pattern)

**What the logical process should do:** Before the ARL begins iterative refinement, assess whether sufficient information exists to operate the loop without returning to the human for missing context. Classify each information need as available, derivable, or missing. Block progression if critical information is missing.

**When it activates:** At loop entry (iteration 0) and at re-entry after a "return to front-loading" decision.

**Why it's needed:** Without information completeness assessment, the loop proceeds with gaps and the agent hallucinate-fills missing information — producing fluent but wrong artifacts. [Q2: SAM ssf:108, :1051; GSD discuss-phase failure mode]

**What success looks like:** The loop begins with all critical prerequisites satisfied. The agent has: the skill's stated purpose, the assessment findings, the plan, and knowledge of downstream references. No runtime escalation for missing context.

**What failure looks like:** The loop begins with missing information and either (a) produces wrong output that must be reverted, or (b) escalates to the human for information that could have been captured upfront.

**Framework patterns that inform this:**

| Pattern | Source | What to Import | What to Adapt |
|---|---|---|---|
| RT-ICA (dynamic prerequisite assessment) | SAM ssf:374-378 | The AVAILABLE/DERIVABLE/MISSING classification model | SAM's RT-ICA is one-shot; the ARL needs re-triggerable RT-ICA within the loop |
| 3-scale front-loading (project/phase/task) | GSD discuss-phase.md:26 | Hierarchical context capture at different granularities | GSD's context is for planning; ARL's context is for skill refinement |
| Fast-path bypass (4 binary criteria) | Octocode research:237-259 | Binary criteria for determining whether deep front-loading is needed | Octocode's criteria are for code review; ARL needs criteria for skill refinement |

**Novel element (not in any framework):** Re-triggerable front-loading within an iterative loop. SAM's front-loading has no "return to discovery" path from execution. The ARL loops — iteration N may discover that iteration 0's assessment was insufficient. The logical process needs a decision point: "Is the current information sufficient to continue, or must we return to the human for additional context?"

---

## R2 — Loop Detection

**What the logical process should do:** Detect when the loop is in a failure state — oscillating (fix A breaks B, fix B breaks A), stalling (same findings recurring), or exceeding resource bounds. Stop or escalate rather than continuing indefinitely.

**When it activates:** At the start of each iteration (before assessment), by comparing current state against prior iteration states.

**Why it's needed:** Without loop detection, the ARL can run indefinitely without converging. The human currently terminates stalled loops by observation — this judgment must be replaced by a detectable condition. [Q3: universal gap — all frameworks lack R2; research doc gap #7]

**What success looks like:** The loop detects its own failure states and either self-corrects (if the failure is recoverable) or escalates with a diagnosis (oscillation detected between states X and Y; stall detected with finding Z recurring for N iterations).

**What failure looks like:** The loop runs for many iterations without converging, consuming resources and producing no net improvement. The human notices the loop is stuck only when they check on it.

**Framework patterns that inform this:**

| Pattern | Source | What to Import | What to Adapt |
|---|---|---|---|
| Task-based thrashing detection | Ralph mod.rs:447, :1917 | Count-based failure escalation (N failures = abandon) | Ralph counts task-level failures; ARL needs finding-level pattern detection |
| Bounded iteration count | GSD plan-phase.md:314-318 | Maximum iteration limit with human decision at boundary (force/guide/abandon) | GSD bounds planning iterations; ARL bounds refinement iterations |
| MaxRuntime cutoff | Ralph mod.rs:429-439 | Time-based termination as safety net | Import directly as a resource bound |

**Novel element (not in any framework):** Finding-level oscillation detection. No framework detects that iteration N's findings resemble iteration N-2's findings. The logical process needs: (1) characterize findings at each iteration, (2) compare across iterations, (3) classify pattern (converging, stalling, oscillating). This requires cross-iteration state — a capability that must be designed from scratch.

**Multi-agent opportunity:** A dedicated loop-monitoring agent could track findings across iterations independently of the working agents, providing an external perspective on loop health. [Q5: multi-agent monitoring pattern]

---

## R3 — Validity Filtering

**What the logical process should do:** Before acting on assessment findings, evaluate their validity. A finding without verifiable evidence (file:line citation, structural observation) is low-confidence. Low-confidence findings consume iteration budget and risk introducing regressions if acted upon.

**When it activates:** After assessment produces findings, before planning acts on them.

**Why it's needed:** The assessor produces findings of varying confidence. Acting on all findings equally wastes iterations on false positives and can introduce regressions. [Q1: BMAD adversarial-review.md:33-35; research doc gap #8]

**What success looks like:** High-confidence findings (with structural evidence) proceed to planning. Low-confidence findings are either verified through a second independent check or deprioritized. False positive rate is measurably lower than unfiltered assessment.

**What failure looks like:** All findings treated equally. False positives trigger changes that must be reverted. Iteration budget consumed by phantom issues.

**Framework patterns that inform this:**

| Pattern | Source | What to Import | What to Adapt |
|---|---|---|---|
| GAN-inspired adversarial validation | Octocode MANIFEST.md:117-125 | Generator/Discriminator tension; cross-model validation | Octocode's pattern validates plans and code; ARL validates assessment findings |
| Forced adversarial review with git truth source | BMAD code-review/instructions.xml:7-14, :26-31 | Cross-referencing claims against external truth source (file system state) | BMAD checks code claims against git; ARL checks findings against skill file structure |
| Two-layer adversarial chain | Ralph ralph.yml:61-63, :100-112 | Reviewer validated by a second reviewer (Handler spot-checks Confessor) | Import the meta-validation pattern; the ARL's "reviewer of the reviewer" |
| Backpressure evidence parsing | Ralph event_parser.rs:320-365 | Machine-verifiable evidence gates (tests pass/fail, structural checks) | Where findings can be verified by structural checks, use machine verification |
| Finding classification (UNCHANGED/UPDATED/INCORRECT) | Octocode PR-reviewer:546-550 | Classify findings before acting; delete INCORRECT findings | Import the classification model for ARL assessment findings |

**Novel element:** Cross-examination between parallel validators. No framework supports the producer disputing a reviewer finding with counter-evidence. The ARL could use two independent assessment agents whose findings are cross-compared — findings confirmed by both are high-confidence; findings from only one are flagged for additional scrutiny. [Q5: cross-examination gap analysis]

---

## R4 — Plan Quality Gates

**What the logical process should do:** Before executing a refinement plan, validate it for internal consistency, proportionality, and alignment with assessment findings. A plan that is inconsistent, addresses the wrong findings, or proposes disproportionate changes should be revised or escalated.

**When it activates:** After planning produces a refinement plan, before implementation begins.

**Why it's needed:** Without plan validation, internally inconsistent or disproportionate plans proceed to implementation, producing changes that must be reverted or that cause more harm than the findings they address. [Research doc gap #9]

**What success looks like:** Plans that pass the quality gate are internally consistent, address the actual findings, and propose changes proportional to finding severity. Failed plans are revised (up to a bounded number of iterations) or escalated.

**What failure looks like:** A plan that rewrites large sections of a skill to address a minor formatting issue proceeds to implementation unchallenged.

**Framework patterns that inform this:**

| Pattern | Source | What to Import | What to Adapt |
|---|---|---|---|
| 7-dimension plan checking with thresholds | GSD gsd-plan-checker.md:55-293 | Multi-dimensional plan evaluation; scope sanity thresholds; adversarial checking instructions | GSD checks execution plans; ARL checks refinement plans. Add proportionality dimension. |
| Adversarial plan verification in GAN flow | Octocode MANIFEST.md:91-94 | Verifier critiques plan for gaps, assumptions, clarity against research evidence | Octocode verifies against research; ARL verifies against assessment findings |
| Bounded revision loop (max 3 iterations) | GSD plan-phase.md:269-312 | Max iteration limit with human decision at boundary | Import directly — prevent unbounded plan revision |

**Parallel opportunity:** Multiple independent agents could check different plan dimensions simultaneously — feasibility, evidence backing, scope proportionality, downstream impact — then synthesize a unified quality verdict. [Q5: parallelizable gate analysis]

---

## R5 — Purpose Anchor

**What the logical process should do:** Record the skill's stated purpose at iteration 0 (before any changes). At each subsequent iteration, compare the skill's current state against the original purpose. Detect when successive changes have cumulatively drifted from the original purpose.

**When it activates:** Captured at iteration 0 (before first assessment). Checked at each iteration boundary (after implementation, before next assessment).

**Why it's needed:** Without a purpose anchor, the loop optimizes for assessment metrics rather than the skill's actual intent. Successive iterations can shift a skill away from its original purpose without any single change being obviously wrong. [Q4: SAM ssf:92, human-out-of-loop:329; research doc gap #1]

**What success looks like:** At iteration N, the skill still serves its original stated purpose. Changes have improved quality without shifting direction. If drift is detected, the loop escalates with specific evidence of what shifted.

**What failure looks like:** After 5 iterations of "improvements," the skill has been optimized for assessment criteria but no longer serves its original audience or use case. The human discovers this only upon manual review.

**Framework patterns that inform this:**

| Pattern | Source | What to Import | What to Adapt |
|---|---|---|---|
| 4-tier deviation rules (acceptable vs unacceptable) | GSD gsd-executor.md:83-139 | Explicit classification of change types: Rules 1-3 (acceptable) vs Rule 4 (unacceptable/architectural) | GSD classifies within single execution; ARL classifies across iterations |
| Objective injection (bookend pattern) | Ralph hatless_ralph.rs:21-23 | Injecting the original objective into every iteration context | Ralph injects but never checks; ARL must inject AND compare |
| Goal-backward verification | GSD gsd-verifier.md:1-14 | "Did the outcome achieve the goal?" verification approach | GSD verifies phase goals; ARL verifies skill purpose preservation |
| Stage 7 Final Verification against Discovery | SAM ssf:757-793 | Comparing final output against original intent capture | SAM does this once at end; ARL must do it at each iteration |

**Novel element:** Cumulative drift detection across iterations. No framework detects gradual semantic drift where individually-acceptable changes cumulatively shift purpose. GSD's deviation rules distinguish acceptable from unacceptable changes within a single execution, but not across a series of iterations. The logical process needs a comparison mechanism that operates on the trajectory of changes, not just the latest delta.

---

## R6 — Content-Loss Detection

**What the logical process should do:** After each implementation step, compare the structural inventory of the skill before changes against the inventory after changes. Detect whether semantic units (sections, headings, behavioral blocks, examples) present before are still present after. Reorganization is acceptable; silent deletion is not.

**When it activates:** After implementation completes, before the next iteration's assessment begins.

**Why it's needed:** During refactoring, content can be silently removed. The current pipeline has no mechanism to detect it. An assessment that finds "all sections well-organized" after implementation may not notice that two sections were deleted. [Q4: universal finding — least-addressed requirement; research doc gap #2]

**What success looks like:** After implementation, a structural comparison confirms all semantic units are present (possibly reorganized). If units are missing, the loop flags the specific missing units and either reverts or escalates.

**What failure looks like:** A refactoring iteration "improves organization" by removing 3 sections deemed "redundant" by the implementing agent. No gate catches the removal. The human discovers the loss only when the skill fails to cover a use case it previously handled.

**Framework patterns that inform this:**

| Pattern | Source | What to Import | Limitation |
|---|---|---|---|
| Content structural inventory comparison (concept) | SAM human-out-of-loop:129 | The conceptual category: before/after structural comparison as machine-verifiable check | Identified as a category but not implemented in any framework |
| AC implementation status checking | BMAD code-review/instructions.xml | Checking implementation against a list of expected items | Checks expected implementations, not whether existing content was preserved |
| Artifact existence/substantive/wired verification | GSD gsd-verifier.md:107 | Three-level existence verification (exists → has substance → is connected) | Checks new artifacts, not preservation of existing content |

**Must be designed from scratch:** Content-loss detection is entirely novel. The logical process needs: (1) capture structural inventory before changes (what sections, headings, examples exist), (2) capture structural inventory after changes, (3) compare inventories to detect missing units, (4) distinguish reorganization (unit moved) from deletion (unit absent). No framework provides any of these capabilities.

---

## R7 — Convergence Tracking

**What the logical process should do:** Track the trajectory of the loop across iterations. Determine whether findings are decreasing (converging), stable (stalling), or alternating (oscillating). Determine when remaining findings are not worth the cost of fixing them (diminishing returns).

**When it activates:** At each iteration boundary, after assessment produces findings for the current iteration.

**Why it's needed:** Without convergence tracking, the loop cannot determine whether it's making progress. It may continue fixing low-value issues indefinitely, or it may oscillate without ever converging. The human currently makes this judgment by observing the loop — this must be replaced by tracked state. [Q3: universal gap; research doc gap #3]

**What success looks like:** The loop tracks finding counts and characteristics across iterations. When findings decrease below a value threshold, the loop terminates with "converged." When findings stall or oscillate, the loop escalates with the specific pattern detected.

**What failure looks like:** The loop runs indefinitely, fixing progressively more trivial issues. Or the loop oscillates between two states, with each "fix" creating a new finding that triggers the previous fix to be undone.

**Framework patterns that inform this:**

| Pattern | Source | What to Import | Limitation |
|---|---|---|---|
| Score delta tracking | GSD gsd-verifier.md, verify-work:480-490 | Tracking previous_score → score across verification runs | Within single phase only; doesn't classify convergence/oscillation/stall |
| Fresh-context-per-iteration (anti-pattern for convergence) | Ralph AGENTS.md:73-85 | Understanding of what is lost when state is not preserved across iterations | Ralph explicitly trades convergence tracking for context freshness |

**Must be designed from scratch:** Cross-iteration convergence tracking is entirely novel. The logical process needs: (1) count and characterize findings at each iteration, (2) store across iterations (requires persistent state), (3) classify trajectory (converging if count decreasing and no oscillation, stalling if count stable, oscillating if pattern detected), (4) apply diminishing returns — when fixing remaining findings costs more than the value they provide, stop.

**Relationship to R2:** R7 (convergence tracking) provides the data that R2 (loop detection) acts on. R7 observes the trajectory; R2 decides what to do about it (continue, escalate, terminate).

---

## R8 — Proportionality Check

**What the logical process should do:** For each finding in the assessment, evaluate whether the proposed fix is proportional to the finding's severity. A low-severity finding (e.g., minor style inconsistency) should not trigger a high-scope change (e.g., rewriting entire section). Flag disproportionate changes for review.

**When it activates:** During plan quality gate (R4), as one dimension of plan evaluation.

**Why it's needed:** Without proportionality checking, the loop acts on all findings equally. Low-severity findings can trigger high-scope changes that introduce risk and consume iteration budget without proportional benefit. [Q3: universal gap; research doc gap #4]

**What success looks like:** The plan's proposed changes are proportional to finding severity. Low-severity findings get low-scope fixes. High-severity findings get whatever scope is necessary. The human sees a proportionality assessment alongside the plan.

**What failure looks like:** A finding that a section heading could be more descriptive triggers a plan to rewrite the entire section's content, organization, and examples.

**Framework patterns that inform this:**

| Pattern | Source | What to Import | Limitation |
|---|---|---|---|
| Scope sanity thresholds | GSD gsd-plan-checker.md:196-208 | Absolute scope limits (tasks/plan, files/plan) as safety bounds | Measures plan SIZE, not proportionality relative to finding severity |
| Change scope classification | BMAD correct-course/instructions.md:161-164 | Classifying changes as Minor/Moderate/Major | Routes to teams, not to proportionality assessment |
| Line count splitting threshold | Octocode plan:149 | Absolute size trigger for decomposition | Measures size, not proportionality |

**Must be designed from scratch:** Proportionality comparison is novel. No framework compares finding severity against proposed change scope. The logical process needs: (1) classify finding severity (impact of not fixing), (2) estimate change scope (how much will be modified), (3) compare — if scope greatly exceeds severity, flag for review or automatically reduce scope.

---

## R9 — Downstream Impact Analysis

**What the logical process should do:** After modifying a skill, identify all components that reference or depend on the modified skill. Verify those references still resolve and contracts still hold. A change that breaks a downstream consumer should be detected before the next iteration.

**When it activates:** After implementation completes, alongside content-loss detection (R6). Both operate as post-implementation verification.

**Why it's needed:** Skills are referenced by other skills, agents, CLAUDE.md, and plugin.json. A refactoring that renames a section, changes a reference path, or removes a capability can break downstream consumers. The assessor's lifecycle audit detects broken references, but only before implementation — not after. [Research doc gap #5]

**What success looks like:** After implementation, all inbound references to the modified skill still resolve. If references are broken, the loop either fixes them (if the fix is proportional) or escalates with the specific broken references.

**What failure looks like:** A refactoring renames a reference file. Three other skills link to the old path. The broken links are not detected until someone invokes those skills.

**Framework patterns that inform this:**

| Pattern | Source | What to Import | Limitation |
|---|---|---|---|
| PR Flow Impact Analysis | Octocode PR-reviewer:491-501 | Tracing downstream effects of changes | Octocode analyzes code changes; ARL analyzes skill file changes |
| Milestone integration checker (cross-phase wiring) | GSD audit-milestone.md:55-71 | Cross-component integration verification | GSD checks at milestone boundary; ARL needs per-iteration checking |
| Redundant convoy observers | Gastown convoy/observer.go:14-67 | Multiple agents monitoring the same lifecycle events | Redundancy for reliability, not reference checking |
| Existing lifecycle audit (bidirectional coherence) | plugin-creator assessor | Reference verification between skills, agents, and plugin.json | Runs before implementation, not after |

**What must be designed:** Post-implementation reference re-verification. The existing assessor lifecycle audits check bidirectional coherence, but they run as part of assessment (before changes), not after implementation. The logical process needs to re-run reference verification after each implementation step. This is an adaptation of an existing capability (lifecycle audit) to a new timing point (post-implementation).

**Parallel opportunity:** Independent agents could verify different reference types simultaneously — markdown links, skill activation references, plugin.json paths, CLAUDE.md references. [Q5: parallel dependency checking]

---

## R10 — Split Justification

**What the logical process should do:** When the refinement loop determines that content belongs in a separate skill, verify that the split produces independently viable units before proceeding. A new skill must be invocable in multiple distinct contexts, not only from its parent. If the new entity is only useful as a child of its parent, it is a reference file, not a skill.

**When it activates:** When the plan proposes extracting content into a new skill (split operation). Before the split is implemented.

**Why it's needed:** Without split justification, the loop can fragment a working skill into pieces that are individually incomplete. The existing `/refactor-skill` splits based on line count and domain boundaries but doesn't verify independent viability. [Research doc gap #6; Q3: universal gap]

**What success looks like:** A proposed split is evaluated against independent viability criteria. If the new skill can be invoked from multiple contexts and serves a standalone purpose, the split proceeds. If the new skill is only meaningful as a subset of its parent, it's classified as a reference file instead.

**What failure looks like:** A 600-line skill is split into three 200-line skills. Two of the three are only ever invoked as references from the original skill — they add navigation complexity without independent value.

**Framework patterns that inform this:**

| Pattern | Source | What to Import | Limitation |
|---|---|---|---|
| Context budget splitting | GSD gsd-plan-checker.md:196-208 | Splitting by resource constraint (plan too large) | Splits by necessity, not by independent viability |
| Size-based splitting | Octocode plan:149 | Triggering split when size exceeds threshold | Measures size, not whether pieces are independently useful |
| Epic/story decomposition | BMAD workflow steps | Domain-based decomposition into units of work | Business domain decomposition, not skill viability |

**Must be designed from scratch:** Independent viability assessment is novel. No framework validates that split products are independently useful. The logical process needs: (1) define "independently viable" for skills (invocable from multiple contexts, serves standalone purpose, understandable without parent), (2) evaluate proposed split against these criteria before implementation, (3) if criteria fail, reclassify the extraction as a reference file rather than a skill.

---

## Cross-Cutting: The ARL Loop Structure

Based on the R1–R10 mapping, the logical process for one ARL iteration follows this sequence of gates:

```
Iteration N begins
  │
  ├─ R1: Information completeness check
  │    (Is context sufficient to proceed? If not, return to human.)
  │
  ├─ R5: Purpose anchor comparison
  │    (Has the skill drifted from its iteration-0 purpose?)
  │
  ├─ ASSESS: Run assessment (existing assessor capability)
  │
  ├─ R3: Validity filtering
  │    (Filter findings — remove low-confidence / false positives)
  │
  ├─ R7 + R2: Convergence and loop state check
  │    (Are findings decreasing? Oscillating? Stalling?)
  │    (If stalling or oscillating → escalate or terminate)
  │    (If diminishing returns → terminate with "converged")
  │
  ├─ PLAN: Generate refinement plan
  │
  ├─ R4 + R8: Plan quality and proportionality gate
  │    (Is the plan internally consistent? Proportional to findings?)
  │
  ├─ R10: Split justification gate (if plan includes splits)
  │    (Are proposed new skills independently viable?)
  │
  ├─ IMPLEMENT: Execute the plan
  │
  ├─ R6: Content-loss detection
  │    (Are all semantic units preserved?)
  │
  ├─ R9: Downstream impact verification
  │    (Do all references still resolve?)
  │
  └─ Loop back to Iteration N+1
```

**Gate timing matters:** R1 and R5 run BEFORE assessment (ensuring the loop is operating in the right context). R3, R7, R2 run AFTER assessment (filtering and evaluating findings). R4, R8, R10 run AFTER planning (validating the plan). R6, R9 run AFTER implementation (verifying changes didn't break anything).

**Exit conditions:** The loop terminates when: (a) R7 indicates convergence (findings below value threshold), (b) R2 detects a failure state (oscillation, stall, resource exhaustion), (c) R1 determines information is insufficient and human input is needed, or (d) maximum iteration count reached.
