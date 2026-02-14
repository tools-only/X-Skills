A systematic multi-agent source code survey protocol for extracting patterns from repositories through structured discussion and cross-examination. Produces evidence-based datasets that inform skills and system design.

## Overview

The expert panel methodology uses multiple AI agents as framework specialists to analyze code repositories through structured questions and cross-examination. The orchestrator coordinates discussion, records evidence trails, and synthesizes findings into two outputs: general theory (reusable patterns) and domain-specific mappings (implementation guidance).

The process is designed for building knowledge bases grounded in primary evidence (source code with file:line citations) rather than secondary sources or training data claims.

## When to Use

Apply this methodology when:

- Building datasets for skills or systems that must be grounded in actual implementations
- Studying multiple frameworks to extract cross-cutting patterns
- Performing comparative analysis of how different codebases solve similar problems
- Extracting reusable design principles from domain-specific implementations
- Answering research questions requiring evidence from multiple independent sources

Do NOT use for:

- Single-codebase documentation (use direct source code analysis)
- Exploratory research without defined questions (too expensive, lacks focus)
- Implementation tasks (this is a research method, not a development workflow)

## Core Protocol

The methodology follows a 4-phase structure. Each phase completes before the next begins.

### Phase 1: Discussion

**Input:** Structured questions anchored to defined requirements or research questions
**Process:** Experts respond from source code, cross-examine each other, orchestrator records
**Output:** Q&A file containing all exchanges with source citations

**Mechanics:**

1. Orchestrator poses a question group to all experts
2. Each expert analyzes their assigned repository and responds with source-cited evidence
3. Experts challenge each other's claims (cross-examination)
4. Orchestrator records: question, responses, challenges, resolutions
5. Repeat for next question group

**Completion criteria:**

- All question groups posed
- Each requirement/research question addressed by minimum 3 experts
- Cross-examination occurred (experts challenged each other)
- Q&A file contains full record

### Phase 2: Requirement Mapping

**Input:** Q&A file from Phase 1, defined requirements
**Process:** Map each requirement to discovered framework evidence
**Output:** Mapping document linking requirements to framework capabilities

**Mechanics:**

For each requirement, document:

- Which frameworks address it? How?
- Which frameworks do NOT address it? What gap does that reveal?
- Where do frameworks disagree on approach? Under what conditions?
- What must be designed from scratch (no framework precedent)?

**Completion criteria:**

- Each requirement has a mapping entry
- Each entry cites specific Phase 1 Q&A exchanges
- No implementation artifacts (schemas, thresholds, code) present

### Phase 3: Synthesis

**Input:** Phase 1 Q&A file, Phase 2 mapping
**Process:** Orchestrator writes synthesis documents drawing only from prior phases
**Output:** General theory document, domain-specific mapping document

**Mechanics:**

1. Write general theory: universal patterns applicable beyond the domain
2. Write domain-specific mapping: which mechanisms from which frameworks inform the target system
3. No new claims introduced at this stage (synthesis only)

**Completion criteria:**

- Both documents exist
- Each section traces to Phase 1 Q&A evidence
- Scope check passes (what/why described, not how to implement)

### Phase 4: Validation

**Input:** Phase 3 synthesis documents, all prior artifacts
**Process:** Stress-test findings for scientific rigor
**Output:** Traceability matrix, limitations section, contribution statement, reproducibility check

**Mechanics:**

1. **Traceability audit:** For every claim, trace backwards through: synthesis → mapping → Q&A → source code. Mark or remove untraced claims.
2. **Research question coverage:** For each research question, document: full answer / partial answer / no answer with section references.
3. **Limitations and threats to validity:** Document single-source risks, inference risks, unresolved tensions, sample boundaries, temporal boundaries.
4. **Contribution statement:** State what is novel — what this process discovered not previously documented.
5. **Reproducibility check:** Verify all artifacts exist for independent replication (repository paths + commits, protocol document, Q&A file, mapping, synthesis).

**Completion criteria:**

- No untraced claims remain (removed or marked unsupported)
- Research question coverage documented
- Limitations section written
- Contribution statement written
- Reproducibility artifacts verified

## Evidence Standards

Evidence hierarchy from strongest to weakest:

| Level     | Type                       | Strength | Examples                                           | Use When                              |
| --------- | -------------------------- | -------- | -------------------------------------------------- | ------------------------------------- |
| Primary   | Source code with file:line | Highest  | "Lines 137-154 of sync.py implement loop detection" | Making definitive claims              |
| Secondary | Documentation, READMEs     | Medium   | "Docs describe RT-ICA as dynamic pattern"         | Describing intent vs implementation   |
| Tertiary  | Cross-framework comparison | Medium   | "3 of 6 frameworks use goal-backward verification" | Identifying design trade-offs         |
| Tertiary  | Experimental observation   | Low      | "Tested 20 tools, accuracy improved 65% to 89%"   | Validating effectiveness empirically  |
| Invalid   | Training data claims       | None     | "Ralph probably uses fuzzy matching"              | Never — must cite actual source       |

**Inadmissible evidence:**

- Claims from model training data without source citation
- Pattern-matched assumptions ("frameworks usually do X")
- Undocumented inferences ("this suggests they intended Y")

**When evidence conflicts:**

- Prefer primary over secondary (code over docs)
- Prefer implemented over documented (working code over design notes)
- Record both positions if unresolved, do not silently pick one side

## Roles

### Orchestrator

**Responsibilities:**

- Pose structured questions to experts
- Record all Q&A exchanges to durable file
- Synthesize findings from evidence (not from own analysis)
- Maintain traceability chain

**Prohibited actions:**

- Answering framework questions itself (experts discover, orchestrator coordinates)
- Delegating synthesis to sub-agents before discussion completes
- Producing implementation artifacts (schemas, thresholds, code)
- Jumping from spawn to file production (skipping discussion phase)

### Experts

**Responsibilities:**

- Discover from source code with file:line citations
- Respond to questions with evidence-based claims
- Challenge other experts' claims (cross-examination)
- Distinguish observable fact from inference

**Prohibited actions:**

- Writing summary files as primary output (discussion is the output)
- Working in isolation without seeing other experts' responses
- Claiming features exist without source code verification
- Fabricating explanations from training data

## Question Design

Effective questions follow these patterns:

### Anchor to Requirements

Each question should advance at least one defined requirement or research question. Questions without this anchor lead to scope creep.

**Example:**

```text
Requirement R3: Validity filtering / false positive handling

Question: How does your framework distinguish real findings from false positives before acting on them?
```

### Use Themed Groups

Group related questions to build context incrementally. Experts can reference prior responses within the group.

**Example theme:** Convergence detection

```text
1. How does your framework detect when iteration is converging vs diverging?
2. What metrics or signals trigger convergence assessment?
3. Under what conditions does your framework declare convergence achieved?
```

### Push with 5 Whys

Surface answers often describe behavior without explaining rationale. Push deeper.

**Example:**

```text
Expert: "We use consensus voting among 3 agents"
Orchestrator: "Why 3 agents instead of 2 or 4?"
Expert: "2 cannot break ties, 4 is more expensive without accuracy gain"
Orchestrator: "Why does tie-breaking matter for your use case?"
Expert: "Binary accept/reject decision must resolve, cannot escalate all ties to human"
```

### Distinguish Fact from Inference

Questions should separate what the code does (observable) from why it does it (inference).

**Example:**

```text
Observable: "What does the code do when zero findings are returned?"
Inference: "Why was this approach chosen over alternative X?"
```

## Cross-Examination Protocol

Cross-examination surfaces contradictions, unimplemented features, and framework-specific assumptions.

### How Cross-Examination Works

After an expert responds, other experts can:

- Challenge factual claims ("I don't see that in the source code at that line")
- Question generalization ("That's true for simple tasks, what about complex tasks?")
- Offer counter-examples ("Our framework does X differently under condition Y")
- Request clarification ("What do you mean by 'soft halt' vs 'hard halt'?")

### Recording Disagreements

When experts disagree and disagreement is not resolved:

- Record both positions
- Note the evidence each cites
- Mark as unresolved tension
- Do NOT silently pick one side

Disagreement is data, not noise. It reveals framework-specific assumptions or contextual factors.

**Example:**

```text
Question: How does your framework handle zero findings?

bmad-expert: "We soft-halt and re-examine (pause, do not block)"
gsd-expert: "We hard-halt (cannot proceed without findings)"

Cross-examination:
bmad-expert: "Why hard-halt instead of re-examining assumptions?"
gsd-expert: "Zero findings means goal already achieved or query was wrong. Both require human decision."

Resolution: Unresolved — depends on whether zero findings is success condition or error condition. Context-dependent.
```

## Context Management

Expert panel sessions consume significant context. These guidelines prevent context exhaustion.

### Incremental Q&A File Writing

The orchestrator must write the Q&A file after each question group, not at the end of the session. The Q&A file is the durable record if context compacts mid-session.

### Phase Boundaries as Breakpoints

If context runs low, Phase 1 can split across sessions. The Q&A file is the handoff artifact — a new session reads it to know what has been discussed and what remains.

### Expert Context Load

Each expert operates in its own context. Give experts focused questions (one requirement or question group at a time) rather than the full question list at once.

### What to Preserve on Compaction

If orchestrator context compacts:

- Q&A file path
- Which questions posed, which requirements mapped
- Which phase the process is in
- Repository paths and expert assignments

## Traceability Chain

Every claim in synthesis must trace backwards through this chain:

```text
Claim in synthesis
  ↓
Section in Phase 2 mapping
  ↓
Q&A exchange in Phase 1
  ↓
Source code citation (file:line)
```

If any link is missing, the claim is unsupported.

**Example of valid traceability:**

```text
Synthesis claim: "Cross-agent consensus reduces false positives more reliably than single-agent confidence scoring"

Phase 2 mapping: "R3 validity filtering — Novel pattern: 3 experts (bmad, gsd, ralph) independently reported multi-agent validation. No single framework implements this, but all three use variations."

Phase 1 Q&A:
  - bmad-expert (lines 245-260 of validator.py): "We run 3 validators in parallel and require 2/3 agreement"
  - gsd-expert (lines 88-92 of consensus.rs): "Findings confirmed by 2+ agents promoted to high-confidence queue"
  - ralph-expert (lines 134-156 of verification.py): "Cross-verification between agents catches hallucinations single agent misses"

Source code: All three citations verified ✓
```

## Known Failure Modes

These patterns cause expert panels to fail. Avoid explicitly.

### Orchestrator Answering Its Own Questions

**Symptom:** Orchestrator analyzes framework repositories and draws conclusions instead of asking experts.

**Why it fails:** Loses the cross-examination benefit. No challenge to orchestrator's interpretation.

**Correct pattern:** Orchestrator poses question, waits for expert responses, synthesizes from their answers.

### Skipping Discussion Phase

**Symptom:** Experts spawned and immediately assigned file-production tasks. No Q&A exchange occurs.

**Why it fails:** Synthesis has no evidence base. Nothing to trace claims back to.

**Correct pattern:** Phase 1 discussion completes before Phase 3 synthesis begins.

### Delegating Synthesis to Agents Without Evidence

**Symptom:** Orchestrator spawns agents to write synthesis documents directly from framework analysis.

**Why it fails:** Agents produce claims without cross-examination or evidence trails. No traceability.

**Correct pattern:** Synthesis draws only from Phase 1 Q&A and Phase 2 mapping.

### File Production Tasks Instead of Questions

**Symptom:** Tasks structured as "write framework-analysis-X.md" instead of "answer question Y."

**Why it fails:** Experts produce isolated summaries. No discussion, no cross-examination.

**Correct pattern:** Tasks are questions. File writing happens after discussion completes.

### Scope Creep into Implementation

**Symptom:** Orchestrator or experts produce YAML schemas, threshold values, file layouts, pseudocode.

**Why it fails:** Implementation details contract solution space. Invented without empirical basis.

**Correct pattern:** Describe what should happen and why (logical process), not how to build it (implementation).

### Not Maintaining Q&A File

**Symptom:** Expert responses scattered across messages. No central record.

**Why it fails:** Loses evidence trail. Cannot trace synthesis claims back to source.

**Correct pattern:** Orchestrator writes Q&A file incrementally after each question group.

## Bias and Limitations

### Selection Bias

The frameworks chosen for analysis do not represent the full landscape. Findings are bounded by the sample.

**Mitigation:** Explicitly document sample (which frameworks, why chosen) and note what populations are excluded.

### Agent Interpretation Bias

Each expert interprets source code through its model's training. Two agents reading the same code may extract different conclusions.

**Mitigation:** Cross-examination partially addresses this. When experts disagree, record both interpretations.

### Absence vs Non-Discovery

When an expert reports "my framework does not implement X," this means the expert did not find it. It does not prove the capability is absent.

**Mitigation:** Distinguish "searched and not found" from "confirmed absent." Note search method and scope.

### Temporal Snapshot

Source code analysis reflects repository state at time of analysis. Frameworks evolve.

**Mitigation:** Note git commit or branch analyzed. Future sessions can check whether findings still hold.

## Reproducibility Requirements

For independent replication, these artifacts must exist:

1. **Repository paths and commits:** Which codebases analyzed, at which git commit
2. **Protocol document:** The question list and process (this methodology)
3. **Q&A file:** The raw evidence from Phase 1
4. **Phase 2 mapping:** The intermediate analysis linking requirements to evidence
5. **Synthesis documents:** The conclusions
6. **Traceability matrix:** The audit trail linking claims to evidence

If any artifact is missing, the synthesis cannot be independently verified. The Q&A file is the most critical — it is the primary data.

## Source Attribution

This methodology was developed through the Autonomous Refinement Loop (ARL) expert panel process (2026-02-12 to 2026-02-13). Six specialist agents analyzed AI development frameworks: BMAD-METHOD, Gastown, Get-Shit-Done, Octocode-MCP, Ralph-Orchestrator, and Stateless Agent Methodology (SAM). The protocol emerged from observing failure modes in prior sessions and designing explicit countermeasures.

**Frameworks analyzed:**

- BMAD-METHOD: `../BMAD-METHOD/`
- Gastown: `../gastown/`
- Get-Shit-Done: `../get-shit-done/`
- Octocode-MCP: `../octocode-mcp/`
- Ralph-Orchestrator: `../ralph-orchestrator/`
- SAM: `./methodology_development/`

**Primary research documents:**

- `plugins/plugin-creator/skills/assessor/references/ARL/autonomous-refinement-loop-research.md`
- `plugins/plugin-creator/skills/assessor/references/ARL/human-out-of-loop-prerequisites.md`
- `plugins/plugin-creator/skills/assessor/references/ARL/ARL-agent-instructions.md`
