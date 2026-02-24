# SDLC Layer Separation Architecture — ARL Integration Suggestions

**Generated**: 2026-02-23  
**Plan**: ARL = Observe → Identify → Probe → Accumulate → Improve  
**Source**: Full read of all ARL items from sdlc-layer-candidates-master.md

---

## Plan Context

- **Observe**: kaizen, session-historian — session recall, transcript analysis
- **Identify**: deviations/hallucinations — fact-check, hallucination-detector, doc-drift
- **Probe**: human for invisible knowledge — (explicit phase; no dedicated tool yet)
- **Accumulate**: local domain knowledge with staleness — research/, KB, context manifest
- **Improve**: SAM + connective tissues — optimize-claude-md, work-backlog-item close, skill/agent patches

---

## 1. agentskill-kaizen (plugin)

**Path**: `plugins/agentskill-kaizen/`

### Integration
- **Observe**: transcript-analyst agent, transcript-analysis skill, kaizen-analysis MCP — session transcript mining
- **Improve**: kaizen-improvement skill, improvement-generator agent — hook generation, agent/skill patches

### Amendment
Add agentskill-kaizen as the **primary Observe→Improve bridge** in the ARL flow. It spans Observe (analysis) and Improve (actionable outputs). The plan should explicitly reference it as the post-hoc implementation of ARL Layer 3 Observation (per plugin-creator/arl).

### Nugget
- **For Probe**: transcript-analysis dimension 9 (Missing Hooks) identifies recurring manual corrections that could be automated — these are candidates for human-probing: "What invisible knowledge makes you correct this manually?"
- **For Accumulate**: kaizen-improvement outputs to `.planning/kaizen/improvements/` — consider routing high-frequency findings into research/ or context manifest as institutional knowledge.

### Extracted Content
- Commands: `/agentskill-kaizen:analyze`, `explore`, `report`, `generate-hooks`
- MCP: kaizen-analysis (process mining, clustering), kaizen-duckdb (persistent DB)
- 10 dimensions: Tool Misuse, Repeated Errors, User Frustration, Missing Tooling, Delegation Patterns, Shortest Path, Red Herrings, System Interruptions, Missing Hooks, DuckDB SQL
- Outputs: hooks, agent patches, skill patches, CLAUDE.md updates, script automation proposals
- Sentiment dashboard: real-time visualization of user message sentiment

---

## 2. transcript-analysis (skill)

**Path**: `plugins/agentskill-kaizen/skills/transcript-analysis/SKILL.md`

### Integration
- **Observe**: Primary skill for session transcript analysis. Place at Observe phase entry point.

### Amendment
Document the JSONL schema and DuckDB query patterns as the canonical **Observe data model**. Plan should specify that Observe phase outputs use `.planning/kaizen/analysis-DATE.md` format for downstream Identify/Improve consumption.

### Nugget
- **For Identify**: Dimension 2 (Repeated Errors) — "Edit-before-Read", "String to replace not found" — these correlate with hallucination-detector's speculation-as-diagnosis; consider cross-referencing.
- **For Probe**: Dimension 3 (User Frustration) — "No,", "Don't", "Why did you" — signals where human knowledge was needed but not captured; feed into Probe design.
- **For Accumulate**: Output format includes "Recommendation type — hook, skill patch, agent prompt fix, CLAUDE.md update" — these map to Accumulate targets (research/, context manifest).

### Extracted Content
- Data: `~/.claude/projects/{project-key}/` — JSONL transcripts, subagent transcripts, tool-results
- Record types: assistant, user, system, progress, file-history-snapshot, summary
- 10 dimensions with extraction methodology (Tool Misuse, Repeated Errors, User Frustration, etc.)
- Process mining: extract_tool_sequences, discover_process_model, check_conformance, find_frequent_patterns, detect_frustration_signals, cluster_sessions
- Output: `.planning/kaizen/analysis-DATE.md` with session ID, severity, evidence, frequency, recommendation type

---

## 3. kaizen-improvement (skill)

**Path**: `plugins/agentskill-kaizen/skills/kaizen-improvement/SKILL.md`

### Integration
- **Improve**: Transform Observe findings into actionable improvements. Place at Improve phase, downstream of Observe.

### Amendment
The plan's "Improve (SAM + connective tissues)" should explicitly include kaizen-improvement as the **transcript-derived improvement pipeline**. It produces delegation prompts for subagent-refactorer, skill-creator, and CLAUDE.md updates — all connective tissue targets.

### Nugget
- **For Probe**: "Improvements are instruction sets for specialist agents, not direct edits" — the delegation protocol could include a Probe step: "Before implementing, ask: does the human have context this finding doesn't capture?"
- **For Accumulate**: Priority scoring includes "Blast radius — project-wide > single-agent > single-session" — project-wide improvements (CLAUDE.md) should update Accumulate artifacts (e.g., research/ patterns).

### Extracted Content
- Prerequisite: Analysis findings in `.planning/kaizen/` from transcript-analysis
- Five improvement types: Hook generation, Agent prompt refinement, Skill patches, CLAUDE.md updates, Script automation
- Delegation protocol: outcome-focused, never prescribe specific code
- Output: `.planning/kaizen/improvements/` (draft) or direct install (hooks only, --install)
- Priority: frequency × impact, automation potential, blast radius, implementation cost

---

## 4. hallucination-detector (plugin)

**Path**: `plugins/hallucination-detector/`

### Integration
- **Identify**: Primary tool for detecting speculation-as-diagnosis, invented causality, pseudo-quantification, completeness overclaims. Place at Identify phase, structural enforcement (Stop hook).

### Amendment
The plan should state that **Identify** includes both (a) skill-based verification (fact-check) and (b) structural enforcement (hallucination-detector Stop hook). The hook blocks completion when triggers are present — architectural constraint, not behavioral instruction.

### Nugget
- **For Probe**: When blocked, Claude must rewrite with evidence-first language. The rewrite often surfaces "I don't know yet" — a natural Probe trigger: "What would you need to know?"
- **For Accumulate**: Triggers (speculation, causality, pseudo-quantification, completeness) could be added to research/ as "hallucination pattern taxonomy" for skill documentation.
- **For Improve**: After 2 blocks in same response cycle, plugin allows completion to prevent infinite loops — document this as a known trade-off in Improve phase.

### Extracted Content
- Stop hook: blocks completion when triggers found in last assistant message
- 5 triggers: Speculation ("I think", "probably"), Causality ("because" without evidence), Pseudo-quantification ("8/10", "70%"), Completeness ("all files checked"), Delegation Micromanagement (overly prescriptive edits)
- Required fix: evidence-first language, cited observations, "I don't know yet"
- Ignores: code blocks, blockquotes, questions
- Trade-off: 2 blocks → allow completion (infinite loop prevention)

---

## 5. plugin-creator/arl (skill)

**Path**: `plugins/plugin-creator/skills/arl/SKILL.md`

### Integration
- **Cross-cutting**: ARL theory reference — informs all phases. Not a phase-specific tool; use when designing or evaluating ARL flows.

### Amendment
The plan should **explicitly reference** the ARL skill as the theoretical foundation. Add a "Theory" subsection that cites: HOOTL, Three Layers (Research Body, Execution Model, Observation), 10 Gates (R1-R10), Scope-Feasibility Matrix, Decision Tree for gate replacement.

### Nugget
- **For Probe**: Layer 2 Execution Model — "Asynchronous feedback queue", "AI user representatives", "Question-to-action-item conversion" — these are Probe design patterns. Plan should add Probe phase using these concepts.
- **For Observe**: "agentskill-kaizen is the current implementation of Layer 3 Observation in post-hoc mode" — plan should state Observe phase = agentskill-kaizen + session-historian.
- **For Accumulate**: Scope-Feasibility Matrix — "scope-classification step must precede any attempt at autonomous operation" — Accumulate should store scope classifications for work items.

### Extracted Content
- HOOTL: human-in-the-loop outcome quality with human-out-of-the-loop execution
- Three Layers: Research Body (R1-R10, 7 principles), Execution Model (pre-discovery, async queue, AI reps, question-to-action), Observation (passive agents, agentskill-kaizen)
- 10 Gates: R1 Information Completeness, R2 Loop Detection, R3 Validity Filtering, R4 Plan Quality, R5 Purpose Anchor, R6 Content-Loss, R7 Convergence, R8 Proportionality, R9 Downstream Impact, R10 Split Justification
- Decision Tree: 4 conditions for replacing human gate (external truth, single dimension, distinguishable success/fail, bounded scope)
- Scope-Feasibility Matrix: High/Medium/Low scope clarity × goal measurability × data enumeration → eliminable?

---

## 6. session-historian (skill)

**Path**: `.claude/skills/session-historian/SKILL.md`

### Integration
- **Observe**: Session recall, transcript search. Place at Observe phase for "I forgot what happened" and context reconstruction.

### Amendment
The plan should distinguish **session-historian** (recall/search) from **transcript-analysis** (mining/patterns). session-historian = user-initiated recall; transcript-analysis = systematic analysis. Both consume `~/.claude/projects/` JSONL.

### Nugget
- **For Probe**: Fidelity rules — "Distinguish absence: 'Not mentioned in transcript' not 'didn't happen'" — when Probe surfaces invisible knowledge, record it so future sessions don't assume absence = didn't happen.
- **For Accumulate**: Summaries cached at `~/.claude/kaizen/session-summaries/` — consider feeding high-value summaries into context manifest or research/ as "session-derived knowledge."
- **For Identify**: session_query.py `search` reads raw JSONL — fact-check and Identify phase can use same index for claim verification against prior session evidence.

### Extracted Content
- Script: `.claude/skills/session-historian/scripts/session_query.py` — list, messages, search, show, index
- Index: `~/.claude/kaizen/session-index.duckdb` (sessions + user_messages tables)
- Summaries: `~/.claude/kaizen/session-summaries/` — AI-generated, mark-summarized
- Fidelity rules: Read before summarizing, verbatim user messages, preserve counts, distinguish absence
- Workflow: list → messages (verbatim) → search (raw) → show + summary template

---

## 7. knowledge-explorer (skill)

**Path**: `.claude/skills/knowledge-explorer/SKILL.md`

### Integration
- **Accumulate**: Manages research/ KB — local domain knowledge with staleness (verified, next_review). Place at Accumulate phase entry point.

### Amendment
The plan should define **Accumulate** as: research/ KB (knowledge-explorer) + context manifest (context-refinement) + Integration Opportunities (research-context-agent). knowledge-explorer is the KB CRUD layer.

### Nugget
- **For Probe**: Staleness (next_review, 6 months) — when KB entry is stale, Probe could ask: "Has anything changed in your experience with this tool?"
- **For Improve**: fetch-github, add, update-append — these operations feed Improve when topic-specialist or research-curator produce verified content.
- **For Identify**: fact-check can verify claims in research/ entries; knowledge-explorer's verified/next_review fields support Identify's evidence discipline.

### Extracted Content
- Script: `research/knowledge-explorer.py` — list, show-template, fetch-github, add, update-append, migrate
- KB root: `research/` with category subdirs
- Frontmatter: name, description, metadata (topic, category, source_url, verified, next_review, tags)
- Staleness: next_review < today or Last Verified > 6 months
- Valid categories: agent-frameworks, developer-tools, mcp-ecosystem, etc. (20+)

---

## 8. refresh-research (skill)

**Path**: `.claude/skills/refresh-research/SKILL.md`

### Integration
- **Accumulate**: Bulk refresh of research entries. Place at Accumulate phase for staleness-driven updates.

### Amendment
The plan should include **refresh-research** as the bulk Accumulate refresh workflow. It runs RT-ICA pre-flight, spawns research-curator agents in waves of 5, updates README, produces summary report.

### Nugget
- **For Probe**: RT-ICA pre-flight checks "mcp__Ref and mcp__exa available" — Probe phase could add: "Human available for invisible-knowledge questions?"
- **For Improve**: Post-actions (lint, commit) — ensure Improve phase's work-backlog-item close and optimize-claude-md don't conflict with refresh-research commits.
- **For Identify**: Staleness detection (past review date, >6 months) — Identify could flag research entries with REFUTED claims from fact-check.

### Extracted Content
- Scope: --all, --stale (default), --category, --dry-run
- Workflow: Inventory → Staleness → RT-ICA → Waves of 5 research-curator → README update → Summary report → Lint + commit
- RT-ICA conditions: mcp__Ref, mcp__exa, gh, network, ./research/ writable
- Outcome categories: Updated, Unchanged, Failed

---

## 9. research-curator (skill)

**Path**: `.claude/skills/research-curator/SKILL.md`

### Integration
- **Accumulate**: Add/maintain research entries. Single URL, batch (--batch), rerun (--rerun), validate (--validate). Place at Accumulate phase as the content-creation layer.

### Amendment
The plan should specify research-curator as the **Accumulate content executor**. It spawns @research-curator agents; refresh-research orchestrates bulk; knowledge-explorer provides CRUD/exploration.

### Nugget
- **For Probe**: Batch mode waves of 5 — when research-curator fails, Probe could ask: "What context from your experience would help verify this?"
- **For Improve**: Post-actions (README, lint, commit, push) — same pattern as work-backlog-item; ensure consistent artifact conventions.
- **For Identify**: Validate mode runs validate_research.py — script detects, agent fixes. fact-check could feed REFUTED claims into --validate workflow.

### Extracted Content
- Modes: Default (single URL), Batch (--batch, waves of 5), Rerun (--rerun), Validate (--validate)
- Agent: @research-curator for content work
- Post-actions: README update, lint, commit, push
- Batch: max 5 concurrent, duplicate detection, progress reporting

---

## 10. research-curator/references (entry-template, validation-rules, batch-mode)

**Path**: `.claude/skills/research-curator/references/entry-template.md`, `validation-rules.md`, `batch-mode.md`

### Integration
- **Accumulate**: Define Accumulate artifact format and validation. Reference from plan's Accumulate schema.

### Amendment
The plan should include the **research entry schema** (entry-template) and **validation rules** (error/warning/info severity) as Accumulate artifact conventions. batch-mode documents wave spawning (max 5) — reuse across ARL phases.

### Nugget
- **For Probe**: Entry template "Relevance to Claude Code Development" — Probe could populate "Invisible knowledge from user experience" subsection.
- **For Improve**: Freshness Tracking (Last Verified, Next Review) — Improve phase (optimize-claude-md, topic-specialist) should update these when modifying research entries.
- **For Identify**: validation-rules — script vs agent responsibility (script detects, agent fixes content) — same pattern as fact-check's claim extraction vs verification.

### Extracted Content
- Entry template: Category selection flowchart, Required Information (Identity, Substance, Relevance), Freshness Tracking
- Validation: section_completeness, header_fields, empty_sections (error); access_dates, freshness_tracking, statistics_currency (warning)
- Batch: waves of 5, sequential, duplicate detection, post-batch README/lint/commit

---

## 11. fact-check (skill)

**Path**: `.claude/skills/fact-check/SKILL.md`

### Integration
- **Identify**: Verify claims against primary sources. Produces VERIFIED/REFUTED/INCONCLUSIVE. Place at Identify phase for backlog items, skill docs, plugin content.

### Amendment
The plan should define **Identify** as: (a) hallucination-detector (structural, Stop hook), (b) fact-check (verification, primary sources), (c) doc-drift-auditor (code vs docs). fact-check is the verification workflow.

### Nugget
- **For Probe**: INCONCLUSIVE verdict — "State what additional step would resolve it" — natural Probe trigger: ask human for invisible knowledge that could resolve.
- **For Accumulate**: VERIFIED/REFUTED claims update backlog; consider routing to research/ when claim is about a tool/library (update KB entry).
- **For Improve**: fact-check spawns @fact-checker in waves of 5 — same pattern as research-curator; topic-specialist consumes fact-check findings for skill updates.

### Extracted Content
- Evidence: WebFetch, WebSearch, CLI output, repo source, MCP — NOT training data
- Claim extraction: backlog item, plugin path, --all-unverified
- Verification: waves of 5 @fact-checker, CoVe (2-3 falsification questions, independent check)
- Verdict: VERIFIED, REFUTED, INCONCLUSIVE with citation
- Post-actions: update backlog, lint, commit

---

## 12. optimize-claude-md (skill)

**Path**: `.claude/skills/optimize-claude-md/SKILL.md`

### Integration
- **Improve**: Optimize AI-facing files (CLAUDE.md, SKILL.md, agents). Place at Improve phase as the "connective tissues" optimizer.

### Amendment
The plan's "Improve (SAM + connective tissues)" should explicitly list optimize-claude-md as the **AI-facing documentation optimizer**. It runs RT-ICA pre-check, 8 optimization principles, CoVe post-check, independent verification.

### Nugget
- **For Probe**: Phase 4 — "If agent signals BLOCKED: Present blocking reason, ask for resolution" — BLOCKED often means missing human context; add Probe step before re-delegation.
- **For Accumulate**: Optimized files are Accumulate targets; ensure optimize-claude-md doesn't remove content that was added from Accumulate (e.g., research-derived patterns).
- **For Identify**: Independent verification (Phase 5) — second agent checks for regressions; aligns with fact-check's CoVe and Identify's evidence discipline.

### Extracted Content
- Phases: Validate → Baseline → Delegate @contextual-ai-documentation-optimizer → Handle BLOCKED/DONE → Independent verification → Measure → Report → Apply on approval
- 8 principles: Positive framing, Motivation, Concrete examples, Front-loaded priorities, Concise language, Explicit format control, Strategic XML tagging, Structural enforcement
- RT-ICA pre-check, CoVe post-check
- Iterative mode for >300 lines: Structural → Content → Polish

---

## 13. work-backlog-item close path (Step 9)

**Path**: `.claude/skills/work-backlog-item/SKILL.md` (Step 9)

### Integration
- **Improve**: Verification agent for acceptance criteria; closes loop on completed items. Place at Improve phase as the **completion verification** step.

### Amendment
The plan should include work-backlog-item close as the **Improve phase completion gate**. It verifies checklist 100%, spawns acceptance-criteria verification agent, writes closing record, closes GitHub issue.

### Nugget
- **For Probe**: Step 9d verification agent prompt — "Does the implementation satisfy the stated goal?" — add: "Were any invisible requirements discovered during implementation?" (feeds Probe/Accumulate).
- **For Accumulate**: Closing record includes "verified by checklist + acceptance criteria check" — context-refinement's "Discovered During Implementation" could be invoked before close to capture Accumulate-worthy findings.
- **For Identify**: Verification agent checks git log, reads changed files — same evidence discipline as fact-check and doc-drift-auditor.

### Extracted Content
- Trigger: $0 = close or resolve
- Close path: 9a Find item → 9c Checklist verification (100%) → 9d Spawn verification agent (PASS/FAIL) → 9e Write closing record
- Verification agent: Read plan, git log, key files; assess goal satisfaction; return PASS/FAIL + evidence
- Resolve path: No verification; reason required
- GitHub: Close issue with comment on success

---

## 14. doc-drift-auditor (agent)

**Path**: `.claude/agents/doc-drift-auditor.md`

### Integration
- **Identify**: Identifies doc vs implementation drift. Place at Identify phase alongside fact-check and hallucination-detector.

### Amendment
The plan should add **doc-drift-auditor** to Identify phase. It complements fact-check (claims) and hallucination-detector (speculation) with **structural drift** — implemented but undocumented, documented but unimplemented, outdated, mismatched details.

### Nugget
- **For Probe**: "Documented but unimplemented" — Probe: "Was this intentionally deferred? What would need to change to implement it?"
- **For Accumulate**: DOCUMENTATION_DRIFT_AUDIT.md findings — feed into context manifest "Discovered During Implementation" or research/ when drift reveals tool/library behavior changes.
- **For Improve**: Drift findings are Improve inputs — update docs to match code, or vice versa; optimize-claude-md could consume drift report for targeted optimization.

### Extracted Content
- Process: Repository discovery → Git timeline → Implementation analysis → Documentation claims → Drift detection → Report
- Categories: Implemented but undocumented, Documented but unimplemented, Outdated, Mismatched details
- Evidence: file:line, commit SHA, quoted claims, code reality
- Output: DOCUMENTATION_DRIFT_AUDIT.md
- Boundaries: Audit only; no modifications; no training-data reliance

---

## 15. code-review (agent)

**Path**: `.claude/agents/code-review.md`

### Integration
- **Identify**: LLM-slop and hallucination detection in code. Place at Identify phase for code artifacts (overlaps with hallucination-detector for narrative, code-review for code).

### Amendment
The plan should include **code-review** in Identify phase for **code-level** deviations: reimplemented scaffolding, junk patterns, placeholders, hallucinated defaults, duplicate env vars. Distinct from hallucination-detector (narrative) and fact-check (claims).

### Nugget
- **For Probe**: "Creating defaults/fallbacks that are entirely hallucinated or imagined" — when code-review flags this, Probe: "What is the correct default for this context?"
- **For Accumulate**: code-review findings (Critical/Warning/Suggestion) — high-frequency patterns could become research/ entries or CLAUDE.md rules.
- **For Improve**: code-review runs during context compaction or pre-commit — ensure Improve phase (work-backlog-item close) invokes code-review before verification agent when code was changed.

### Extracted Content
- Focus: "Some or all of the code was generated by an LLM"
- Patterns: Reimplemented scaffolding, junk patterns, placeholders, hallucinated defaults, duplicate env vars, indentation/JSON/YAML issues
- Severity: Critical (security, correctness), Warning (reliability, performance), Suggestion (alternatives, docs)
- Process: Get changes → Understand patterns → Focus areas → Review against standards
- Boundaries: Focus on what matters; respect existing choices; be specific

---

## 16. context-refinement (agent)

**Path**: `.claude/agents/context-refinement.md`

### Integration
- **Accumulate**: Updates context manifest with discoveries from work session. Place at Accumulate phase as the **institutional knowledge** capture agent.

### Amendment
The plan should define context-refinement as the **Accumulate session-capture** agent. It reads transcript, identifies drift/discoveries, appends "Discovered During Implementation" to context manifest. Reduces staleness of context for future sessions.

### Nugget
- **For Probe**: "Wrong assumptions in original context" — when context-refinement finds this, Probe: "What assumption should have been in the original context? How would we capture it next time?"
- **For Observe**: context-refinement reads transcript from `sessions/transcripts/context-refinement/` — ensure this aligns with session-historian and transcript-analysis data locations.
- **For Improve**: "Guardian of institutional knowledge" — Improve phase (work-backlog-item close) should invoke context-refinement before closing to capture last-session discoveries.

### Extracted Content
- Trigger: End of work session
- Process: Read transcript → Analyze for drift/discoveries → Decision (no update vs update) → Append "Discovered During Implementation"
- Qualifies: Undocumented interactions, incorrect assumptions, hidden side effects, complex error cases
- Doesn't qualify: Minor typos, implied things, temporary workarounds
- Output: "No context updates needed" or "Context manifest updated with X discoveries"

---

## 17. logging (agent)

**Path**: `.claude/agents/logging.md`

### Integration
- **Observe**: Consolidates work session output. Place at Observe phase as the **task state maintainer** for future session consumption.

### Amendment
The plan should include **logging** in Observe phase. It ensures task file reflects present state for kaizen/session-historian consumption. "Maintains clean task state for future sessions" — enables Observe to read accurate work logs.

### Nugget
- **For Accumulate**: logging's Work Log (Completed, Decisions, Discovered, Next Steps) — context-refinement could read this to avoid duplicating "Discovered" in context manifest.
- **For Improve**: logging runs during context compaction or task completion — coordinate with work-backlog-item close so logging runs before close verification.
- **For Probe**: "Decisions" and "Discovered" sections — when these are sparse, Probe: "Were there decisions or discoveries not captured?"

### Extracted Content
- Trigger: Context compaction or task completion
- Responsibilities: Read file, read transcript, assess cleanup, remove irrelevant, update existing, add new, chronological order
- Format: Completed, Decisions, Discovered, Next Steps
- Rules: Cleanup first, chronological integrity, consolidation, clarity
- Boundaries: Only edit task file; never touch sessions/state/, current-task.json

---

## 18. topic-specialist (agent)

**Path**: `.claude/agents/topic-specialist.md`

### Integration
- **Improve**: Researches primary sources, can update or create skills with verified findings. Place at Improve phase as the **verified-knowledge-to-skill** bridge.

### Amendment
The plan should include **topic-specialist** in Improve phase. It feeds SAM + connective tissues with primary-source-verified content. OUTPUT: "answer + update skill" or "answer + create skill".

### Nugget
- **For Probe**: "If a primary source is unavailable: state 'Unable to verify from primary source'" — Probe: "Do you have experience with this that could substitute?"
- **For Accumulate**: topic-specialist populates skills from verified findings — skills are Accumulate targets; ensure research/ entries are updated when topic-specialist adds tool/library knowledge.
- **For Identify**: topic-specialist uses fact-check, find-cause, research-curator skills — CoVe, evidence discipline; aligns with Identify's verification requirements.

### Extracted Content
- Invocation: TOPIC, SKILLS, QUESTION, CONDITIONS, OUTPUT (answer only | answer + update skill | answer + create skill)
- Research: GitHub source, README, docs, issues — primary sources only
- CoVe: 2-3 falsification questions, cross-check
- Update skill: Append with citations, do not remove existing
- Create skill: Invoke skill-creator, add-doc-updater, populate with verified findings
- Boundaries: No training-data-only; no commit; no plugin/agent creation

---

## 19. research-context-agent (agent)

**Path**: `.claude/agents/research-context-agent.md`

### Integration
- **Accumulate**: Cross-references research with skills, agents, hooks, commands. Discovers integration opportunities. Place at Accumulate phase as the **connective tissue** between research KB and capabilities.

### Amendment
The plan should include **research-context-agent** in Accumulate phase. It appends "Integration Opportunities" to research files — Enhances Existing, New Skill Candidates, New MCP Candidates, Cross-References. This is the Accumulate→Improve bridge (identifies what to improve).

### Nugget
- **For Probe**: "New skill candidate" — when research-context-agent proposes one, Probe: "Is there existing institutional knowledge that would inform this skill?"
- **For Improve**: Integration Opportunities table — direct input to Improve phase (topic-specialist, optimize-claude-md, kaizen-improvement).
- **For Identify**: research-context-agent validates claims against primary sources (WebSearch/WebFetch) — aligns with fact-check evidence discipline.

### Extracted Content
- Process: Absorb (extract from research) → Search & Match (5 dimensions) → Append (Integration Opportunities)
- Dimensions: Enhance skills, enhance agents, enhance hooks, enhance commands, new skill candidate, new MCP candidate
- Output: Enhances Existing table, New Skill Candidates, New MCP Server Candidates, Cross-References
- Rules: Concrete over vague, skip empty sections, no false positives, preserve content, idempotent

---

## Summary: ARL Phase Mapping

| Phase | Skills | Agents | Plugins |
|-------|--------|--------|---------|
| **Observe** | session-historian, transcript-analysis | logging, transcript-analyst | agentskill-kaizen |
| **Identify** | fact-check | doc-drift-auditor, code-review, fact-checker | hallucination-detector |
| **Probe** | (none — design needed) | (none) | — |
| **Accumulate** | knowledge-explorer, refresh-research, research-curator | context-refinement, research-context-agent, research-curator | — |
| **Improve** | optimize-claude-md, kaizen-improvement, work-backlog-item (close) | topic-specialist, improvement-generator | agentskill-kaizen |
| **Cross-cutting** | plugin-creator/arl | — | plugin-creator |

---

## Recommended Plan Amendments

1. **Add Probe phase design**: Use ARL Layer 2 (async feedback queue, AI user representatives, question-to-action-item) and transcript-analysis dimension 3 (User Frustration) to design human-probing workflow.
2. **Unify wave spawning**: Document "max 5 concurrent" as ARL-wide pattern (fact-check, research-curator, refresh-research, groom-backlog-item).
3. **Define Accumulate schema**: research/ entry template + context manifest "Discovered During Implementation" format.
4. **Define Improve completion sequence**: logging → context-refinement → work-backlog-item close (with verification agent) → code-review (if code changed).
5. **Reference plugin-creator/arl**: Add "ARL Theory" subsection citing HOOTL, 10 Gates, Scope-Feasibility Matrix, Layer 3 Observation (agentskill-kaizen).
