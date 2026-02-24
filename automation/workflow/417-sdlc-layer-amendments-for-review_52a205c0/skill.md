# SDLC Layer Plan — Amendments for Review

**Purpose**: Key amendments from the four integration suggestion reports, extracted for review before merging into the plan. Approve what you want merged; reject or defer the rest.

---

## Layer 0 Amendments

| # | Amendment | Source | Impact |
|---|-----------|--------|--------|
| L0-1 | **development-harness as canonical Layer 0 hub** — Add explicit statement: "The development-harness plugin is the reference implementation of the SAM 7-stage pipeline. It owns: process orchestration, ARL human touchpoint decisions, artifact management, and fallback behavior." | L0 #1 | Structural |
| L0-2 | **Loop limits** — Add as Layer 0 constants: NEEDS_WORK loop limit = 3 iterations per task; NOT_CERTIFIED loop limit = 2 iterations before escalation | L0 #2 | Deliverable |
| L0-3 | **Escalation format** — When escalating, present: (1) what stage, (2) what triggered it, (3) what the agent knows, (4) what the agent does not know, (5) decision options. Never present vague "please review" requests | L0 #3 | Deliverable |
| L0-4 | **Seven artifact types** — List DISCOVERY, PLAN, CONTEXT, TASK, EXECUTION, REVIEW, VERIFICATION with required sections per type | L0 #4 | Deliverable |
| L0-5 | **planner-rt-ica vs rt-ica** — Planner RT-ICA allows planning under uncertainty (APPROVED-WITH-GAPS); execution RT-ICA blocks on missing. Any task under APPROVED-WITH-GAPS MUST pass rt-ica before execution | L0 #8 | Clarification |
| L0-6 | **verification-gate 4 checkpoints** — Run before S5 Execution self-verification: (1) Hypothesis stated, (2) Hypothesis verified, (3) Hypothesis-action alignment, (4) Pattern-matching detection. Activate before any Bash/Write/Edit/NotebookEdit | L0 #11 | Deliverable |
| L0-7 | **orchestrator-discipline anti-patterns** — Investigation Escalation (3+ Read/Grep/Bash without Edit/Write/Task), Agent Output Polling (TaskOutput block=false). "Orchestrators delegate; agents implement. No exemption categories." | L0 #12 | Deliverable |
| L0-8 | **Fact-check before RT-ICA** — REFUTED → MISSING. If fact-check returns REFUTED on a claim, RT-ICA should treat it as MISSING | L0 | Clarification |
| L0-9 | **Remove verification-gate "5% cost for 95% reliability"** — Marked REFUTED (unsubstantiated). Use qualitative phrasing instead | L0 #11 | Amendment |

---

## Layer 1 Amendments

| # | Amendment | Source | Impact |
|---|-----------|--------|--------|
| L1-1 | **language-manifest-template as canonical starting point** — Require all Layer 1 language plugins to produce a manifest conforming to this template before composing with the harness. Distinguish: template = quick-start skeleton; schema = validation rules | L1 #1 | Deliverable |
| L1-2 | **Linting Discovery Protocol** — Extract as Layer 1 standard: "Language plugins MUST implement a discovery sequence before executing quality gates." (git hook → CI → fallback) | L1 #2 | Deliverable |
| L1-3 | **Workflow pattern taxonomy** — Each language plugin should document: TDD-equivalent, Feature Addition, Code Review, Refactoring, Debugging — with agent chains and quality gates per pattern | L1 #3 | Deliverable |
| L1-4 | **typecheck: (none)** — For non-typed languages (e.g. Bash), allow `typecheck: (none)` in quality gates | L1 | Schema |
| L1-5 | **Harness role ↔ agent archetype mapping table** — Add explicit mapping from abstract roles (architect, test-designer, etc.) to agent archetypes | L1 | Deliverable |
| L1-6 | **CoVe bypass anti-pattern** — "Orchestrators pass paths and outcomes; agents discover and verify." Do not pre-gather data for agents | L1 #3 | Clarification |

---

## Layer 2 Amendments

| # | Amendment | Source | Impact |
|---|-----------|--------|--------|
| L2-1 | **Stack profiles list** — Add to plan: documentation-authoring, gitlab, github, mcp_server_development, structured_data, c_cpp_formatting, plugin_development, backlog_management, milestone_workflow, changelog_generation, daily_releases, browser_automation | L2 | Deliverable |
| L2-2 | **STATUS block as shared output contract** — DONE, BLOCKED, FAILED. Consider standardizing across stack profiles (STATUS, SUMMARY, ARTIFACTS, VALIDATION, NOTES) | L2 #1 | Schema |
| L2-3 | **Stack profile extends for composition** — e.g. daily-releases extends changelog_generation. Support inheritance/composition in stack profile schema | L2 | Schema |
| L2-4 | **Workflow definition schema** — From the-rewrite-room frontmatter: required [workflow, version, output_contract], optional [canonical_agent, canonical_skill, canonical_path], validation_gates [HARD_STOP, SOFT_STOP] | L2 #2 | Schema |

---

## ARL Amendments

| # | Amendment | Source | Impact |
|---|-----------|--------|--------|
| A-1 | **Probe phase design** — No dedicated Probe tool yet. Design should use ARL Layer 2 Execution Model concepts: async feedback queue, AI user representatives, question-to-action-item conversion. Use transcript-analysis dimension 3 (User Frustration) as Probe trigger | ARL #5, #2 | Design |
| A-2 | **Phase mapping** — Observe (kaizen, session-historian, logging) → Identify (hallucination-detector, fact-check, doc-drift-auditor, code-review) → Accumulate (knowledge-explorer, refresh-research, research-curator, context-refinement, research-context-agent) → Improve (kaizen-improvement, optimize-claude-md, work-backlog-item close, topic-specialist) | ARL | Structural |
| A-3 | **Identify = structural + verification** — Identify includes (a) hallucination-detector (structural, Stop hook), (b) fact-check (VERIFIED/REFUTED/INCONCLUSIVE), (c) doc-drift-auditor (code vs docs) | ARL #4, #11 | Clarification |
| A-4 | **Accumulate definition** — research/ KB (knowledge-explorer) + context manifest (context-refinement) + Integration Opportunities (research-context-agent). knowledge-explorer = KB CRUD layer | ARL #7 | Clarification |
| A-5 | **Probe trigger from Dimension 3** — User Frustration signals ("No,", "Don't", "Why did you") = where human knowledge was needed but not captured. Feed into Probe design | ARL #2 | Design |
| A-6 | **Probe trigger from Dimension 9** — Missing Hooks = recurring manual corrections. Candidates for human-probing: "What invisible knowledge makes you correct this manually?" | ARL #1 | Design |
| A-7 | **ARL theory subsection** — Add explicit "Theory" subsection citing: HOOTL, Three Layers, 10 Gates, Scope-Feasibility Matrix, Decision Tree. Note: 4 build-from-scratch requirements (R6, R7, R8, R10) emerge from iterative refinement — invisible to single-pass | ARL #5 | Deliverable |
| A-8 | **session-historian vs transcript-analysis** — session-historian = user-initiated recall/search; transcript-analysis = systematic mining/patterns. Both consume ~/.claude/projects/ JSONL | ARL #6 | Clarification |

---

## Summary of Amendments by Type

| Type | Count | Layers |
|------|-------|--------|
| Structural | 4 | L0, L1, ARL |
| Deliverable | 14 | L0, L1, L2, ARL |
| Schema | 4 | L1, L2 |
| Clarification | 6 | L0, L1, ARL |
| Design | 4 | ARL |
| Amendment (remove/reframe) | 1 | L0 |

---

## Recommended Merge Order

1. **Clarifications** (A-8, L0-5, L1-6, A-3, A-4) — Low risk, improve accuracy
2. **Structural** (L0-1, A-2, A-7) — Define hubs and phase mapping
3. **Deliverables** (L0-2, L0-3, L0-4, L0-6, L0-7, L1-1, L1-2, L1-3, L2-1) — Add to plan sections
4. **Schema** (L1-4, L1-5, L2-2, L2-3, L2-4) — Add to Layer 2 schema design
5. **Design** (A-1, A-5, A-6) — Add to Probe phase design (Open Question 3)
6. **Amendment** (L0-9) — Remove refuted claim

---

## How to Use

1. Review each amendment. Mark ✅ Approve, ❌ Reject, or ⏸ Defer.
2. Return edits to the plan with approved amendments merged.
3. Deferred items can be added to Open Questions or a future backlog.
