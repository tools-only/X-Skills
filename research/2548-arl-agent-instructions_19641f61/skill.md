# ARL Agent Instructions

**Purpose:** Context and instructions for agents working on the Autonomous Refinement Loop (ARL) for the plugin-creator's assessor skill. Use this document to understand the project, avoid past mistakes, and execute correctly.

---

## 0. Scope Boundary — Logical Process Only

**This pass produces the logical process design. It does NOT produce implementation.**

The goal is to define:

- **What** must happen at each stage of the loop
- **When** each check, gate, or decision occurs relative to other stages
- **Why** each element is needed — what failure mode it prevents, what it enables
- **Decision trees** — under what conditions does the loop take path A vs path B?
- **Expected outcomes** — what does success look like at each gate? What does failure look like?
- **Goals and standards** — what quality bar must be met? How is "done" defined?
- **Frameworks and methodologies** — which patterns from which frameworks inform each element?

**Referencing existing implementations is allowed.** Describing how Framework X achieves something is evidence — it informs the logical process design. The boundary is between "Ralph detects thrashing by comparing consecutive task descriptions" (observation, referenceable) and "the ARL will detect thrashing by comparing consecutive task descriptions using fuzzy matching at 90% threshold" (premature commitment to a specific approach). Record what exists. Do not prescribe what the ARL will do at the implementation level.

The goal is NOT to define:

- **How the ARL will implement** any of this (no YAML schemas, no file paths, no pseudocode, no algorithm specifications for the ARL)
- **Threshold values for the ARL** (no `0.90`, no `max_lines: 50`, no iteration counts)
- **Data structures for the ARL** (no `.arl/state/` directories, no artifact token formats)
- **Tool choices for the ARL** (no `tiktoken`, no `rapidfuzz`, no embedding model selection)

**Why this boundary matters:** Implementation details invented without evidence contract the solution space. A YAML schema written now locks in assumptions about storage format. A threshold value written now has no empirical basis. The logical process must be clean and clear first — what should happen, when, and why. Implementation is a separate pass that happens later, informed by the logical process.

**Enforcement:** If the orchestrator or any expert begins producing implementation artifacts (schemas, code, file layouts, config formats), stop and redirect to the logical question: "What should happen here and why?"

**Prior session violation:** The previous session's synthesis documents contain implementation details — YAML schemas, proportionality matrices with specific line-count limits, convergence algorithm pseudocode, file path conventions. None of that came from expert discussion or source code evidence. Treat those sections as scope violations to be replaced with logical process descriptions.

---

## 1. What Is ARL?

The **Autonomous Refinement Loop** is a loop that runs without human babysitting:

```
Assess → Plan → Implement → Review → Repeat until converged
```

**Goal:** Replace human judgment at each gate with machine-verifiable conditions so the loop can run autonomously (or with human-at-end-only).

**Primary research files (read these first):**

- `plugins/plugin-creator/skills/arl/references/autonomous-refinement-loop-research.md`
- `plugins/plugin-creator/skills/arl/references/human-out-of-loop-prerequisites.md`

---

## 1a. Research Questions

The expert panel process is structured around these research questions. All discussion, mapping, and synthesis should advance answers to them.

**Primary research question:**

> Under what conditions can human judgment at iterative refinement loop gates be replaced by machine-verifiable conditions, and what are the failure modes when those conditions are insufficient?

**Secondary research questions:**

1. What patterns for autonomous loop control exist across the surveyed frameworks, and which are implemented (not merely documented)?
2. For each ARL requirement (R1–R10), what is the minimum set of logical conditions that must hold for the gate to operate without human intervention?
3. Where do the surveyed frameworks disagree on approach, and what contextual factors (task complexity, bounded vs unbounded scope, risk level) explain those disagreements?
4. What capabilities does no surveyed framework address, and what does that absence reveal about the state of autonomous iterative refinement?
5. What general principles emerge that apply beyond the ARL — to any autonomous development loop?

**How research questions guide the process:** Every Phase 1 discussion question (Section 6) should advance at least one research question. Every Phase 2 mapping entry should cite which research question it informs. Every Phase 3 synthesis section should trace to a research question. If a section cannot be traced to a research question, it is out of scope.

---

## 1b. Methodology

This section describes the expert panel process as a research method — how evidence is gathered, validated, and synthesized. This is not a description of ARL's internal process; it is a description of how we study ARL.

**Method: Multi-agent source code survey with cross-examination**

Six specialist agents each analyze one framework's source code repository. The orchestrator poses structured questions (Section 6) anchored to defined requirements (R1–R10). Experts respond with source-cited evidence. Experts cross-examine each other's claims. The orchestrator records all exchanges in a Q&A file.

**Why this method:** The frameworks exist as code repositories, not as published papers. Traditional literature review cannot access them. Agent-based source code analysis allows systematic survey of implementations. Cross-examination between experts surfaces contradictions, unimplemented features, and framework-specific assumptions that a single analyst would miss.

**Evidence standards:**

- **Primary evidence:** Source code with file:line citations. This is the strongest form — the code either does or does not do what is claimed.
- **Secondary evidence:** Framework documentation, README files, design notes. These describe intent, which may or may not match implementation. Always verify against primary evidence.
- **Tertiary evidence:** Cross-framework comparison. When multiple frameworks solve the same problem differently, the disagreement itself is evidence of a design trade-off worth analyzing.
- **Inadmissible:** Claims from training data, pattern-matched assumptions, undocumented inferences. If it cannot be cited, it cannot be used.

**Bias and limitations of the method:**

- **Selection bias:** The six frameworks were chosen for proximity and relevance to Claude Code plugin development. They do not represent the full landscape of autonomous development tools. Findings are bounded by this sample.
- **Agent interpretation bias:** Each expert agent interprets source code through its own model's training. Two agents reading the same code may extract different conclusions. Cross-examination partially mitigates this, but does not eliminate it.
- **Absence vs non-discovery:** When an expert reports "my framework does not implement X," this means the expert did not find it. It does not prove the capability is absent. The Q&A file should distinguish "searched and not found" from "confirmed absent."
- **Temporal snapshot:** Source code analysis reflects the state of each repository at the time of analysis. Frameworks evolve. Findings should note the commit or branch analyzed where possible.

**Reproducibility requirements:**

For this process to be reproducible, a future session (or independent researcher) needs:
1. The same six repository paths, at the same commits
2. This instructions document (the research protocol)
3. The Q&A file from Phase 1 (the raw evidence)
4. The Phase 2 mapping (the intermediate analysis)

If any of these are missing, the synthesis cannot be independently verified. The Q&A file is the most critical — it is the primary data.

---

## 2. Original Vision (Two-Part Outcome)

When standing up framework expert teams or doing ARL design work, the desired outcome is:

1. **General theory** — A "general theory of relativity" for autonomous development: reusable best practices, principles, and patterns that apply across scenarios and tiers of complexity. Usable as a template for _any_ autonomous system in the future.

2. **ARL-specific synthesis** — Concrete import mappings: which mechanisms from which frameworks (BMAD, Gastown, GSD, Octocode, Ralph, SAM) directly support the plugin-creator's ARL. What to import, what to adapt, what must be built from scratch.

**Knowledge requirement:** Build a web of knowledge grounded in **real data and source code evidence**, not speculation.

---

## 2a. Prior Work State

Previous sessions produced synthesis documents and intermediate findings. This section records what was discovered — not what should be built.

**Git state:** The synthesis documents were committed in prior sessions but later deleted from the working tree. They exist in git history but are not on disk. The content in those documents mixed logical process observations with implementation artifacts (YAML schemas, threshold values, file path conventions). The logical observations are valuable; the implementation artifacts violated scope.

**Key corrections from prior expert panel (verified against source code):**

| Claim (Prior Session) | Correction | Source |
| ---------------------- | ---------- | ------ |
| Ralph uses 90% fuzzy similarity for loop detection | Ralph's loop detection exists only as a Python prototype concept, not in the production Rust codebase. Only task-based thrashing detection is implemented. | ralph-expert source code analysis |
| SAM's RT-ICA is a checklist-based gate | RT-ICA is a dynamic information-completeness assessment pattern — it adapts to the goal, not a fixed checklist | SAM methodology docs |
| BMAD halts hard on zero findings | BMAD's zero-findings response is a soft halt (pause and re-examine) not a hard block | bmad-expert source code analysis |
| GSD verifier uses line-by-line diff | GSD verification is goal-backward (did the outcome occur?) not line-level comparison | gsd-expert source code analysis |
| Gastown tracks convergence via finding count | Gastown does not implement convergence tracking — it uses phase-based progression | gastown-expert source code analysis |

**Novel patterns discovered (not in any single framework):**

1. **Cross-agent validity consensus** — Multiple independent agents evaluating the same finding reduces false positives more reliably than single-agent confidence scoring
2. **Goal-backward verification** — Checking "did the expected outcome occur?" rather than "were all steps executed?" catches process-correct-but-outcome-wrong failures
3. **Tiered autonomy** — Simple/bounded tasks run fully autonomous; complex/unbounded tasks escalate. The distinction is bounded by whether the definition of done is machine-checkable
4. **Structural diff for content loss** — Detecting whether sections, headings, or semantic units disappeared during editing, rather than line-level diff
5. **Purpose anchor at iteration zero** — Recording the stated intent before any changes begin, then checking drift against that anchor at each iteration

**What these findings mean for the expert panel:** The corrections table shows where prior sessions accepted framework claims without source-code verification. The novel patterns show cross-framework insights that no single expert would have surfaced alone. Both should inform the discussion questions — experts should be asked to verify or challenge these findings.

---

## 3. Framework Expert Team

| Expert | Assigned Repository | GitHub URL |
| --- | --- | --- |
| bmad-expert | ../BMAD-METHOD/ | <https://github.com/bmad-code-org/BMAD-METHOD> |
| gastown-expert | ../gastown/ | <https://github.com/steveyegge/gastown> |
| gsd-expert | ../get-shit-done/ | <https://github.com/glittercowboy/get-shit-done> |
| octocode-expert | ../octocode-mcp/ | <https://github.com/bgauryy/octocode-mcp> |
| ralph-expert | ../ralph-orchestrator/ | <https://github.com/mikeyobrien/ralph-orchestrator> |
| sam-expert | ./methodology_development/ | (in-repo) |

See [expert-repos.md](./expert-repos.md) for the full manifest with URL sources.

Experts discover and report what their framework contributes — do not prescribe answers. The orchestrator asks; experts respond from their own analysis of the code.

---

## 4. Anti-Patterns (What NOT to Do)

From prior sessions — avoid these explicitly:

1. **Experts writing summary files** — The user wanted collaborative discussion, not experts producing framework-analysis-\*.md files. Do not have remote experts write files as their primary output.

2. **Experts working in isolation** — Experts must see each other's answers and be able to cross-examine. Discussion should be visible to all participants.

3. **Verifying summaries instead of discussing ARL** — The focal point should be "how would your framework add/change/enhance the ARL?" not "verify and correct the framework analysis files."

4. **Importing unimplemented mechanisms** — Always trace claims to source code. Documentation and design notes can describe features that were never implemented.

---

## 4a. Orchestrator Anti-Patterns

From prior sessions — the orchestrator itself caused these failures:

1. **Answering framework questions itself** — The orchestrator analyzed framework repositories and drew its own conclusions instead of asking experts. Experts exist to discover; the orchestrator exists to coordinate.

2. **Delegating synthesis to sub-agents without discussion** — The orchestrator spawned agents to write synthesis documents directly, skipping the expert discussion phase entirely. Synthesis should follow discussion, not replace it.

3. **Jumping from spawn to verification** — Experts were spawned and immediately assigned file-production tasks. The discussion phase (Section 6 questions) was skipped. The sequence should be: spawn → discuss → synthesize.

4. **Creating tasks for file production instead of questions** — Tasks were structured as "write framework-analysis-X.md" instead of "answer R3: How does your framework handle validity filtering?" File production is an output of the process, not a task for experts.

5. **Scope creep into implementation** — The orchestrator and agents produced YAML schemas, threshold values, file path conventions, and pseudocode. None of this was requested. The scope is logical process design (Section 0).

6. **Losing the Q&A file** — No Q&A file was maintained. Expert responses were scattered across agent messages with no central record. The orchestrator must maintain the Q&A file as discussion proceeds.

---

## 5. arl-expert-panel Command Specification

When creating or improving `.claude/commands/arl-expert-panel.md`, the command must:

1. **Bootstrap a team** of framework experts assigned to their respective repositories.

2. **Orchestrator drives discussion** — The orchestrator asks questions and records Q&A to a file. Remote experts do **not** write files; they respond to questions only.

3. **Prerequisites** — All participants (orchestrator + experts) read the two ARL research docs first so discussion is anchored to the right context.

4. **Cross-examination** — Experts can see each other's answers and ask each other about their areas. The goal is collaborative, critical discourse, not patch-style contributions.

5. **Evidence-based responses** — Experts must respond with validated, unguessed information. Source citations and file:line evidence where applicable.

6. **Discussion anchors** — Anchor questions to ARL requirements (R1–R10) from `autonomous-refinement-loop-research.md` to avoid open-ended meandering.

7. **Tone** — Experts are critical, considerate of desired outcomes, and focused on advancing the ARL's context harness and process for achieving autonomy and high-quality, vetted, reliable, predictable output.

---

## 5a. Agent Teams: Capabilities the Orchestrator Must Use

When running the arl-expert-panel, the lead agent must leverage agent teams as follows. These are implementation requirements, not optional.

**Spawn prompt is the only context teammates receive.** Teammates do not inherit the lead's conversation history. Each spawn prompt must include:

- Paths to the two ARL research docs (so they read them)
- The expert's assigned repository path
- Role: discover and report from source code; respond to questions; do not write files; message other teammates to cross-examine
- That other teammates exist and can be messaged (see `~/.claude/teams/{team-name}/config.json` `members` array for discovery)

**Delegate mode** — Operate as coordination-only. Spawn, message, manage tasks, synthesize. Do not implement or answer framework questions yourself — that is the experts' role. If the system supports delegate mode, use it to enforce this; otherwise self-restrict.

**Task list** — Create tasks the experts will claim. Structure around R1–R10 (e.g. "Answer R3: How does your framework handle validity filtering?") or by phase (e.g. "Initial response on convergence detection" then "Cross-examine bmad-expert's answer"). Aim for 5–6 tasks per teammate so work stays scoped and the lead can reassign if blocked.

**Messaging** — Teammates have `message` (1:1) and `broadcast` (all). Instruct them to use these for cross-examination: after answering, message or broadcast to challenge other experts' claims. The lead receives all messages automatically. Synthesize that stream into the Q&A file — questions asked, to whom, answers, challenges, resolutions.

**Members discovery** — Teammates can read `~/.claude/teams/{team-name}/config.json` to get the `members` array (agent IDs, names). Use this when instructing experts to message specific teammates by name.

---

## 5b. Execution Phases

The expert panel runs in three phases. Each phase completes before the next begins.

**Phase 1: Discussion (the core of the process)**

The orchestrator poses the Section 6 questions to the team. Experts respond from their source code analysis. Experts challenge each other's answers. The orchestrator records all Q&A to the Q&A file.

This phase is the majority of the work. It produces the raw material — verified observations, cross-framework tensions, areas of agreement, areas where frameworks diverge, and areas where no framework addresses the problem.

**Phase 1 completion criteria:** All Section 6 question groups have been posed and discussed. Each R1–R10 requirement has been addressed by at least 3 experts. Cross-examination has occurred — experts have challenged each other's claims. The Q&A file contains the full record.

**Phase 2: Mapping R1–R10 to framework evidence**

Using the Q&A from Phase 1, the orchestrator maps each ARL requirement (R1–R10) to what was discovered:
- Which frameworks address this requirement? How?
- Which frameworks do NOT address it? What gap does that reveal?
- Where do frameworks disagree on approach? What conditions explain the disagreement?
- What must be designed from scratch (no framework covers it)?

This is still analysis, not implementation. The output describes what the logical process should do at each requirement, informed by framework evidence.

**Phase 2 completion criteria:** Each R1–R10 has a mapping entry. Each entry cites specific framework evidence from Phase 1 Q&A. No entry contains implementation artifacts (schemas, code, thresholds).

**Phase 3: Synthesis writing**

The orchestrator writes the two synthesis documents (Section 8) and finalizes the Q&A file. Synthesis draws only from Phase 1 and Phase 2 outputs. No new claims are introduced at this stage.

**Phase 3 completion criteria:** Both synthesis documents exist. Each section traces to Phase 1 Q&A evidence. The scope check passes — every section describes what and why, not how.

**Phase 4: Validation and rigor review**

After synthesis is written, this phase stress-tests the work for scientific rigor. It does not add new findings — it challenges what was produced.

**4a. Traceability audit** — For every claim in both synthesis documents, trace backwards: synthesis section → Phase 2 mapping entry → Phase 1 Q&A exchange → source code citation. If any link in that chain is missing, the claim is unsupported. Mark it as such or remove it.

**4b. Evaluation against research questions** — Review each research question (Section 1a). For each one: does the synthesis provide an evidence-based answer? A partial answer with acknowledged gaps? Or no answer at all? Document the coverage explicitly. Unanswered research questions are not failures — they are findings about the limits of the current evidence.

**4c. Limitations and threats to validity** — Document what could be wrong with the findings:
- Which conclusions depend on a single expert's analysis? (Single-source risk)
- Which conclusions extrapolate beyond what the source code shows? (Inference risk)
- Where did experts disagree and the disagreement was not resolved? (Unresolved tension)
- What would change if a seventh framework were added? (Sample boundary)
- What would change if the same frameworks were analyzed at a different commit? (Temporal boundary)

**4d. Contribution statement** — State explicitly what is novel in the findings. What did this process discover that was not known before? The novel patterns (Section 2a) are candidates, but Phase 4 should validate whether they held up under Phase 1 cross-examination or were revised.

**4e. Reproducibility check** — Verify that all artifacts needed for independent replication exist:
- Repository paths and commits analyzed
- This instructions document (the protocol)
- The Q&A file (the raw data)
- The Phase 2 mapping (the intermediate analysis)
- The synthesis documents (the conclusions)

If any artifact is missing, note it. If the Q&A file is incomplete (e.g., context ran out mid-discussion), note which questions were and were not covered.

**Phase 4 completion criteria:** Traceability audit complete with no untraced claims remaining (removed or marked unsupported). Research question coverage documented. Limitations section written. Contribution statement written. Reproducibility artifacts verified.

---

## 5c. Context Window Management

Expert panel sessions consume significant context. These guidelines prevent context exhaustion before the process completes.

**Orchestrator context load:** The orchestrator receives all expert messages automatically. With 6 experts and cross-examination, message volume is high. The orchestrator should write the Q&A file incrementally (after each question group, not at the end) so that the Q&A file serves as the durable record if context compacts.

**Expert context load:** Each expert operates in its own context. Experts should be given focused questions (one R-requirement or one question group at a time) rather than the full Section 6 list at once.

**Phase boundaries as natural breakpoints:** If context runs low, Phase 1 can be split across sessions. The Q&A file is the handoff artifact — a new session reads the Q&A file to know what has been discussed and what remains.

**What to preserve on context compaction:** The Q&A file path, which questions have been posed, which R-requirements have been mapped, and which phase the process is in. The two ARL research doc paths. The team name and expert assignments.

---

## 6. Discussion Questions for the Expert Panel

The orchestrator poses these questions to the team. Experts discuss, pull them apart, challenge each other, and apply the 5 whys to drill into root causes. Do not accept surface answers — push until the underlying mechanics are clear.

### Human Gating and Judgment

- When working with AI, why does a human stop the AI from doing an action? (5 whys: keep asking why until you reach root causes.)
- What does the human need to know to be sure about what the AI is doing, and what it will do next?
- At what points is judgment needed? Which of those could be eliminated, front-loaded, or reduced?
- What processes need more judgment than others? How do we know?
- How do systems and processes get evaluated as needing closer human attention? What criteria or heuristics do frameworks use?

### Interaction Points and Front-Loading

- What are the interaction points in these frameworks? Where does the human touch the system?
- How does each framework assess what is needed by the human at each of those points?
- How do they front-load human-based work so it can be done all at once — allowing the agent to continue without the human until the task is done?
- What must be captured upfront for that to work? What fails when it isn't captured?

### Framework Design Rationale and Conditions

- Why has your framework chosen that way of doing it? What problem was it solving?
- Is that way the only way, or is it conditional? Under what conditions does the framework do X vs Y?
- How are those conditions detected? What triggers the framework to take one path vs another?
- What makes a project, feature, or task hit these complexity or bounded/unbounded checkpoints? How does the framework decide that something is "simple" vs "complex," or "bounded" vs "unbounded"?

### Alignment, Drift, and Definition of Done

- How is alignment between human intent and agent behavior tracked?
- How is drift detected? What constitutes drift vs acceptable evolution?
- How is the definition of done validated? Who validates it? What makes it machine-checkable vs judgment-dependent?
- When definition of done is ambiguous, how do frameworks handle it?

### Agent Teams and Orchestration

- Given tools like agent teams (separate instances, messaging, shared tasks, lead coordination), how could this be utilized within a system to assist with orchestration and human-gating?
- Where could parallel teammates replace sequential human checks?
- Where could cross-examination between teammates surface issues before they reach a human?
- How could a lead agent synthesize escalation so the human gets compressed context instead of raw logs?

---

## 7. ARL Requirements (R1–R10)

Use these as discussion anchors when querying experts. The two ARL research docs define the problems; experts discover how their framework addresses each.

- **R1** — RT-ICA / front-loading as gate pattern
- **R2** — Loop detection (output-similarity or task-based thrashing)
- **R3** — Validity filtering / false positive handling
- **R4** — Plan quality gates
- **R5** — Purpose anchor (record intent at iteration 0, check drift)
- **R6** — Content-loss detection (before/after diff)
- **R7** — Convergence tracking (finding count per iteration)
- **R8** — Proportionality check (severity vs blast radius)
- **R9** — Downstream impact analysis
- **R10** — Split justification (when to extract a new skill)

---

## 8. Synthesis Outputs (When Running the Process)

When the expert panel process runs successfully, it should produce:

1. **synthesis-general-theory.md** — Universal patterns applicable to any autonomous development system. The expert panel discovers and synthesizes these; structure is not predetermined. Contains: principles, decision trees, when/why rationale, expected outcomes, quality standards. Does NOT contain: schemas, pseudocode, file layouts, threshold values.

2. **synthesis-arl-applicable.md** — ARL-specific mapping of R1–R10 to framework mechanisms: what to import, what to adapt, what must be built from scratch. For each requirement: what the logical process should do, when it activates, why it's needed, what success/failure looks like, and which framework patterns inform it. Does NOT contain: implementation specifications, data structures, algorithm details, config formats.

3. **Q&A file** — The durable record of the discussion. Orchestrator maintains this; experts do not write it. Updated incrementally after each question group (not deferred to end of session).

   **Q&A file structure:**

   For each question posed, record:
   - The question (verbatim)
   - Which experts responded
   - Each expert's response (attributed, with source citations where given)
   - Challenges raised by other experts
   - Resolution or unresolved disagreement
   - Which R-requirements the discussion touches

   The Q&A file is the handoff artifact between sessions. A new session reads it to know what has been discussed and what remains. It is also the evidence base for Phase 2 mapping and Phase 3 synthesis — if a claim appears in synthesis but not in the Q&A file, it has no evidence trail.

4. **Traceability matrix** — A table mapping each synthesis claim to its Phase 2 entry, Phase 1 Q&A exchange, and source code citation. Produced during Phase 4a. This is the audit trail — it allows any claim to be independently verified by following the chain backwards.

5. **Limitations and threats to validity** — Produced during Phase 4c. Documents single-source risks, inference risks, unresolved tensions, sample boundaries, and temporal boundaries. This section makes the work honest about what it does and does not establish.

6. **Research question coverage** — Produced during Phase 4b. For each research question (Section 1a), states whether the synthesis provides a full answer, partial answer, or no answer, with pointers to the relevant synthesis sections.

7. **Contribution statement** — Produced during Phase 4d. States what is novel — what this process discovered that was not previously documented in any single framework.

**Scope check before writing synthesis:** For each section, ask: "Am I describing what should happen and why, or am I describing how to build it?" If the answer is "how to build it," stop and reframe as the logical process.

---

## 9. Infrastructure Gaps (ARL Today)

These gaps describe what the ARL currently lacks. The "What's Missing" column describes the logical capability gap — what should happen but doesn't. The "Prior Panel Finding" column notes what the expert panel discovered about this gap, if anything.

| Gap | What's Missing | Prior Panel Finding |
| --- | -------------- | ------------------- |
| Purpose anchor | Nothing records skill intent at iteration 0 and checks it hasn't drifted | No framework implements this directly. SAM's RT-ICA captures intent at planning time but does not re-check it at each iteration. This is a build-from-scratch capability. |
| Content-loss detection | No before/after structural diff to detect whether semantic units (sections, headings, behavioral blocks) disappeared during editing | GSD's verifier checks goal achievement, not structural preservation. The structural diff concept was identified as a novel cross-framework pattern (Section 2a). |
| Convergence tracking | No mechanism to track whether finding counts are decreasing across iterations | Gastown does not track convergence — prior session incorrectly attributed this. No framework implements iteration-aware convergence tracking. Build-from-scratch. |
| Proportionality check | No evaluation of whether a fix's blast radius is proportional to the severity of the finding it addresses | BMAD's risk assessment is project-level, not finding-level. No framework operates at finding-level proportionality. The logical process needs to define when a fix is "too large" relative to what it corrects. |
| Downstream impact | No post-implementation re-verification of inbound references — after changing a skill, nothing checks whether other skills that reference it still work | Octocode's dependency graph tracks what references what, but does not re-verify after changes. The logical process needs to define when downstream verification activates and what "still works" means. |
| Split justification | No check that new skills created by splitting have 3+ independent trigger scenarios | No framework addresses this. The assessor skill currently checks skill size but not whether a proposed split produces viable independent skills. |
| Loop detection | No detection of whether the loop is producing diminishing returns or oscillating between states | Ralph's task-based thrashing detection is implemented; output-similarity detection is not. The logical process needs to define what constitutes "oscillating" vs "converging slowly." |
| Validity filtering | No mechanism to distinguish real findings from false positives before acting on them | Cross-agent consensus (Section 2a novel pattern #1) was identified as the most promising approach. The logical process needs to define what constitutes agreement and what happens on disagreement. |
| Plan quality gates | No validation that a plan is sound before executing it | BMAD and GSD both have plan review stages. The logical process needs to define what "sound" means for an ARL iteration plan — distinct from project-level plan quality. |

---

## 10. Verification Notes

- **RT-ICA** (Reverse Thinking – Information Completeness Assessment): Use before planning. Confirm prerequisites (framework repos, ARL doc paths, discussion structure) before executing.

- **Source-code verification**: When a framework is claimed to implement X, trace to actual code. Documentation and design notes can describe unimplemented features.

- **Traceability standard**: Every claim in synthesis must trace to Q&A evidence. Every Q&A response must trace to source code. A claim without this chain is unsupported — mark it or remove it.

- **Absence vs non-discovery**: "Expert did not find X" is not the same as "X does not exist." Record what was searched and what was not found. Do not conclude absence from a failed search.

- **Disagreement is data**: When two experts contradict each other on the same topic, the disagreement is itself a finding. Record both positions, the evidence each cites, and whether resolution was reached. Do not silently pick one side.

- **Commit pinning**: Where possible, note the git commit or branch of each framework repository at the time of analysis. This allows future sessions to check whether findings still hold against the current state of the code.

- **Research question alignment**: If a discussion thread, mapping entry, or synthesis section cannot be linked to at least one research question (Section 1a), it is out of scope for this pass. It may be valuable — but it belongs in a future pass, not this one.
