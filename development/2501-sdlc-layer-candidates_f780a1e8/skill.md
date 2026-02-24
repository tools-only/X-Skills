# SDLC Layer Separation Architecture — Skill Candidates

**Generated**: 2026-02-23  
**Source**: Full read of all skills under `.claude/skills/` and their `references/` subdirectories  
**Plan Context**: Layer 0 (SDLC-Agnostic), Layer 1 (Language-Specific), Layer 2 (Stack/Goal-Specific), ARL Meta-Layer

---

## Layer 0 (SDLC-Agnostic)

### work-backlog-item
- **Type**: skill
- **Path**: `.claude/skills/work-backlog-item/SKILL.md`
- **Proposed Layer**: L0
- **Relevance**: Core bridge between BACKLOG.md and SAM planning pipeline; defines human touchpoints (interactive browser, AskUserQuestion), artifact conventions (Plan field, BACKLOG.md structure), and task format for feature requests.
- **Key content**: Step-by-step workflow (find item → groom → RT-ICA gate → compose feature request → invoke SAM → update BACKLOG with Plan reference); `--auto` mode for agent-only execution; close/resolve paths with checklist verification and acceptance-criteria agent spawn; GitHub Issue sync as optional integration.

### work-backlog-item/references/sam-definition.md
- **Type**: reference
- **Path**: `.claude/skills/work-backlog-item/references/sam-definition.md`
- **Proposed Layer**: L0
- **Relevance**: Canonical SAM definition for the repo — stateless agents, externalized memory, RT-ICA gate, artifact tokens, verification at boundaries, human touchpoint model.
- **Key content**: Core principles (stateless agents, message passing, verification at boundaries, RT-ICA gate, semantic artifact tokens); 7-stage pipeline (Discovery → Planning → Context Integration → Task Decomposition → Execution → Forensic Review → Final Verification); artifact flow; escalation rules (3 NEEDS_WORK loops, 2 NOT_CERTIFIED loops).

### work-backlog-item/references/github-integration.md
- **Type**: reference
- **Path**: `.claude/skills/work-backlog-item/references/github-integration.md`
- **Proposed Layer**: L0 (optional integration)
- **Relevance**: Documents BACKLOG.md ↔ GitHub Issue field mapping and human touchpoints for issue creation, milestone assignment, closure.
- **Key content**: Step 2.5 (issue sync), Step 2.5a (create issue), Step 2.7 (in-progress label), Step 9 extension (close issue); setup-github command; story-format body template.

### work-backlog-item/references/validation-plan.md
- **Type**: reference
- **Path**: `.claude/skills/work-backlog-item/references/validation-plan.md`
- **Proposed Layer**: L0
- **Relevance**: Verification protocol for GitHub integration — V1–V6 commands and full integration test sequence.
- **Key content**: V1 (labels), V2 (issue creation), V3 (milestone), V4 (in-progress label), V5 (closure), V6 (BACKLOG.md consistency); full integration test sequence.

### groom-backlog-item
- **Type**: skill
- **Path**: `.claude/skills/groom-backlog-item/SKILL.md`
- **Proposed Layer**: L0
- **Relevance**: Orchestrates backlog grooming with fact-check → RT-ICA → groomer agents; produces context manifests and grooming report. SDLC-agnostic discovery workflow.
- **Key content**: Fact-check before RT-ICA (REFUTED → MISSING); RT-ICA per item; spawn @backlog-item-groomer with RT-ICA + fact-check context; grooming report with Fact-Check Results, RT-ICA Results, Cross-Item Findings.

### rt-ica
- **Type**: skill
- **Path**: `.claude/skills/rt-ica/SKILL.md`
- **Proposed Layer**: L0
- **Relevance**: Mandatory pre-planning checkpoint; blocks planning until prerequisites verified. Core SDLC-agnostic gate.
- **Key content**: Goal reconstruction → reverse prerequisite enumeration → availability verification (AVAILABLE/DERIVABLE/MISSING) → completeness decision (APPROVED/BLOCKED); condition categories (functional, non-functional, interfaces, environment, data, access, operational, delivery, verification, risks); integration with CoVe-style planning; question templates for missing inputs.

### subagent-contract
- **Type**: skill
- **Path**: `.claude/skills/subagent-contract/SKILL.md`
- **Proposed Layer**: L0
- **Relevance**: Global contract for specialist subagents — role boundaries, scope discipline, DONE/BLOCKED status signaling.
- **Key content**: Identity constraints (specialist only, no scope expansion, no invention); behavioral rules (follow SOP, minimal scope, report commands, prefer reversible); DONE/BLOCKED signaling format; quality verification checklist; scope discipline; supervisor protocol; anti-patterns (scope creep, assumption making, silent blocking).

### fact-check
- **Type**: skill
- **Path**: `.claude/skills/fact-check/SKILL.md`
- **Proposed Layer**: L0
- **Relevance**: Primary-source claim verification; evidence discipline; training data recall explicitly rejected.
- **Key content**: Evidence rules (WebFetch, WebSearch, CLI output, repo source, MCP tools valid; training data recall invalid); claim extraction and classification; wave spawning (5 concurrent); verdict format (VERIFIED/REFUTED/INCONCLUSIVE); CoVe requirement (initial lookup → verification questions → independent check → final verdict); report format with citations.

### find-cause
- **Type**: skill
- **Path**: `.claude/skills/find-cause/SKILL.md`
- **Proposed Layer**: L0
- **Relevance**: Evidence-chain protocol for investigations; transforms vague requests into reproducible-proof investigations.
- **Key content**: Evidence chain rules (command output, file content, direct observation valid; docs intent, training recall, inference from absence invalid); Step 0 (capability discovery); Step 1 (disambiguate + success criteria); Step 1.5 (prerequisite check, reproduction safety, bound vs unbound constraints); Step 2 (reproduce and observe); Step 3 (read source); Step 4 (build evidence chain); Step 5 (present findings); prohibited behaviors.

### verify
- **Type**: skill
- **Path**: `.claude/skills/verify/SKILL.md`
- **Proposed Layer**: L0
- **Relevance**: Rigorous self-assessment checklist before marking task complete; prevents premature completion claims.
- **Key content**: Task type & strategy; "WORKS" check (executable vs static); "FIXED" check (reproduction, resolution); quality gates; honesty check; golden rule (demonstrate with evidence); verification summary format.

### delegate
- **Type**: skill
- **Path**: `.claude/skills/delegate/SKILL.md`
- **Proposed Layer**: L0
- **Relevance**: Quick delegation template for sub-agent prompts; WHERE-WHAT-WHY framework.
- **Key content**: Template (OBSERVATIONS, DEFINITION OF SUCCESS, CONTEXT, AVAILABLE RESOURCES, YOUR TASK); delegation rules (Observations + Success Criteria + Resources - Assumptions - Micromanagement); no HOW, constraints OK; checklist.

### cove-prompt-design
- **Type**: skill
- **Path**: `.claude/skills/cove-prompt-design/SKILL.md`
- **Proposed Layer**: L0
- **Relevance**: Chain of Verification for prompt design; reduces hallucinations via generation → verification → revision.
- **Key content**: Three phases (initial answer, verification questions, final validated answer); when to use (factual accuracy, costly hallucinations, multiple facts); verification question design (independent, falsifiable); CoVe vs Chain of Thought; failure modes; practical template.

### scientific-thinking
- **Type**: skill
- **Path**: `.claude/skills/scientific-thinking/SKILL.md`
- **Proposed Layer**: L0
- **Relevance**: Hypothesis-driven reasoning for complex problems; systematic investigation.
- **Key content**: Observation (factual only) → Hypothesis (H₀, Hₐ) → Prediction → Experiment Plan (Tree-of-Thought) → Execution Control (confounds) → Execute and Record → Conclusion; when to use (bug, architecture, refactor, strange behavior); when NOT to use (trivial tasks).

### commit-staged
- **Type**: skill
- **Path**: `.claude/skills/commit-staged/SKILL.md`
- **Proposed Layer**: L0
- **Relevance**: Conventional commits format; artifact convention for git history.
- **Key content**: Type (feat, fix, docs, style, refactor, test, chore); scope rules (WHERE in codebase, not change type); scope anti-patterns; breaking change footer; template workflow.

---

## Layer 1 (Language-Specific)

### agent-creator
- **Type**: skill
- **Path**: `.claude/skills/agent-creator/SKILL.md`
- **Proposed Layer**: L1
- **Relevance**: Abstract roles (Researcher, Planner, Coder, Creator, Tester, Reviewer, etc.); agent schema as language manifest for Claude Code agents.
- **Key content**: Phase 1–7 workflow (Discovery, Requirements, Template Selection, Adaptation, Creation, Validation, Placement); model selection guide; permission mode guide; tool access patterns; description writing guide; role-based contract archetypes; standard vs role-based template decision.

### agent-creator/references/agent-schema.md
- **Type**: reference
- **Path**: `.claude/skills/agent-creator/references/agent-schema.md`
- **Proposed Layer**: L1
- **Relevance**: Complete frontmatter specification; language manifest for agent definitions.
- **Key content**: Required fields (name, description); optional (model, tools, disallowedTools, permissionMode, skills, hooks, maxTurns, mcpServers, memory, background, isolation, color); validation rules; YAML multiline bug.

### agent-creator/references/agent-templates.md
- **Type**: reference
- **Path**: `.claude/skills/agent-creator/references/agent-templates.md`
- **Proposed Layer**: L1
- **Relevance**: Role archetypes as abstract roles; standard vs role-based template methodology.
- **Key content**: Standard templates (user-facing, flexibility); role-based contract archetypes (orchestrated, DONE/BLOCKED); base skeleton; Researcher, Planner/Architect, Coder, Creator, Tester, Reviewer, DevOps/SRE, Auditor, Context Gatherer, Optimizer, Domain Expert; supervisor co-prompt templates; best practices.

### agent-creator/references/agent-examples.md
- **Type**: reference
- **Path**: `.claude/skills/agent-creator/references/agent-examples.md`
- **Proposed Layer**: L1
- **Relevance**: Real-world agent implementations; reference examples for role patterns.
- **Key content**: Plugin assessor, code review, context optimizer, doc drift auditor, skill refactorer, context gathering; role-based examples (Coder Next.js+Supabase, Coder Python TUI, Formatter with hooks); pattern summary; anti-patterns.

### skill-research-process
- **Type**: skill
- **Path**: `.claude/skills/skill-research-process/SKILL.md`
- **Proposed Layer**: L1
- **Relevance**: Systematic process for building skills; quality gates; anti-hallucination checkpoint.
- **Key content**: Stage 1 (Initialize, categorization agent, TODO checklist); Gate 1 (category verification); Stage 2 (parallel category research); Gate 2 (anti-hallucination checkpoint, citation required); Stage 3 (integration); Gate 3 (final validation); MCP tool selection; citation format; agent team alternative for Stage 2.

### skill-research-process/references/agent-prompts.md
- **Type**: reference
- **Path**: `.claude/skills/skill-research-process/references/agent-prompts.md`
- **Proposed Layer**: L1
- **Relevance**: Agent prompt templates for categorization and research; abstract role prompts.
- **Key content**: (Referenced but not read in full — categorization agent, research agent templates.)

### skill-research-process/references/mcp-tools.md
- **Type**: reference
- **Path**: `.claude/skills/skill-research-process/references/mcp-tools.md`
- **Proposed Layer**: L1
- **Relevance**: Tool fidelity guide; project detection for MCP availability.
- **Key content**: (Referenced — WebFetch low, mcp__exa medium, mcp__Ref high.)

### skill-research-process/references/gaps-analysis.md
- **Type**: reference
- **Path**: `.claude/skills/skill-research-process/references/gaps-analysis.md`
- **Proposed Layer**: L1
- **Relevance**: Known limitations; improvement opportunities for skill creation.
- **Key content**: (Referenced — gaps and improvement opportunities.)

### optimize-claude-md
- **Type**: skill
- **Path**: `.claude/skills/optimize-claude-md/SKILL.md`
- **Proposed Layer**: L1
- **Relevance**: Optimizes AI-facing files; 8 optimization principles; quality gates for SKILL.md completeness.
- **Key content**: Phase 1–8 (validate, baseline, delegate to optimizer, handle response, independent verification, measure output, report, apply on approval); 8 principles (positive framing, motivation, concrete examples, front-loaded priorities, concise language, explicit format control, XML tagging, structural enforcement); RT-ICA pre-check, CoVe post-check; iterative mode for large targets; scope expansion (skill dir, plugin dir).

---

## Layer 2 (Stack/Goal-Specific)

### gh
- **Type**: skill
- **Path**: `.claude/skills/gh/SKILL.md`
- **Proposed Layer**: L2
- **Relevance**: GitHub CLI setup and usage; toolchain config for proxy environments; stack-specific (GitHub).
- **Key content**: Installation (setup_gh.py); -R flag for proxy environments; common commands (PR, issues, labels, projects, milestones, workflows, releases); automation script (github_project_setup.py); output formatting.

### gh/references/milestones.md
- **Type**: reference
- **Path**: `.claude/skills/gh/references/milestones.md`
- **Proposed Layer**: L2
- **Relevance**: Milestone CRUD; naming conventions; stack-specific.
- **Key content**: gh api REST; PyGithub; @octokit/rest; automation script; naming conventions.

### gh/references/labels.md
- **Type**: reference
- **Path**: `.claude/skills/gh/references/labels.md`
- **Proposed Layer**: L2
- **Relevance**: Label taxonomy (priority, type, status); BACKLOG.md mapping.
- **Key content**: Priority (p0–idea), type (feature, bug, refactor, docs, chore), status (in-progress, blocked, needs-grooming, needs-review); gh CLI, PyGithub, Octokit; bulk setup.

### gh/references/projects-v2.md
- **Type**: reference
- **Path**: `.claude/skills/gh/references/projects-v2.md`
- **Proposed Layer**: L2
- **Relevance**: GitHub Projects V2; custom fields; GraphQL; stack-specific.
- **Key content**: gh project CLI; GraphQL for field editing; Status (Backlog, Grooming, In Progress, Review, Done); standard project structure.

### gh/references/issue-stories.md
- **Type**: reference
- **Path**: `.claude/skills/gh/references/issue-stories.md`
- **Proposed Layer**: L2
- **Relevance**: Issue as story format; lifecycle; BACKLOG.md ↔ Issue mapping.
- **Key content**: Title convention (conventional commits); body template (Story, Description, Acceptance Criteria, Context); lifecycle (needs-grooming → in-progress → closed); field mapping.

### create-backlog-item
- **Type**: skill
- **Path**: `.claude/skills/create-backlog-item/SKILL.md`
- **Proposed Layer**: L2 (goal-specific: backlog management)
- **Relevance**: Creates BACKLOG.md items; optional GitHub Issue for P0/P1; stack-specific (claude_skills repo).
- **Key content**: Guided, quick, --auto modes; required fields; duplicate detection; compose item block; write to BACKLOG.md; GitHub Issue creation (P0/P1); gh CLI usage.

### create-milestone
- **Type**: skill
- **Path**: `.claude/skills/create-milestone/SKILL.md`
- **Proposed Layer**: L2
- **Relevance**: Creates GitHub milestone; stack-specific (Jamie-BitFlight/claude_skills).
- **Key content**: Guided/quick modes; duplicate check; create via Python script; next steps (group-items-to-milestone, start-milestone).

### create-merge-request-changelog
- **Type**: skill
- **Path**: `.claude/skills/create-merge-request-changelog/SKILL.md`
- **Proposed Layer**: L2
- **Relevance**: MR/PR description generation; domain-based change categorization; works with any git repo.
- **Key content**: analyze_git_changes.py; AI analysis prompts; format_mr_description.py; change categories (bug fixes, enhancements, tech debt, docs, testing, build/CI, non-functional); breaking changes detection; conventional commits integration; GitLab/GitHub CLI integration.

### daily-releases
- **Type**: skill
- **Path**: `.claude/skills/daily-releases/SKILL.md`
- **Proposed Layer**: L2
- **Relevance**: GitHub Releases with AI changelogs per day; reuses create-merge-request-changelog pipeline.
- **Key content**: list_daily_ranges.py; analyze_git_changes per day; AI analysis; format_mr_description; publish_daily_release.py; idempotent (skips up-to-date).

### group-items-to-milestone
- **Type**: skill
- **Path**: `.claude/skills/group-items-to-milestone/SKILL.md`
- **Proposed Layer**: L2
- **Relevance**: Assigns BACKLOG items to milestone; creates missing GitHub Issues; updates Project V2.
- **Key content**: Resolve milestone; load backlog; present selection; create issues; assign existing; update Project V2; report.

### start-milestone
- **Type**: skill
- **Path**: `.claude/skills/start-milestone/SKILL.md`
- **Proposed Layer**: L2
- **Relevance**: Bulk-transitions status labels; updates Project V2 to In Progress.
- **Key content**: Resolve milestone; list issues; confirm; ensure label; bulk label transition; update Project V2.

### complete-milestone
- **Type**: skill
- **Path**: `.claude/skills/complete-milestone/SKILL.md`
- **Proposed Layer**: L2
- **Relevance**: Closes milestone; handles open issues (carry forward, remove, close); updates Project V2 to Done.
- **Key content**: Audit; handle open issues (A/B/C/D options); close milestone; update Project V2; completion report.

### agent-browser
- **Type**: skill
- **Path**: `.claude/skills/agent-browser/SKILL.md`
- **Proposed Layer**: L2
- **Relevance**: Browser automation CLI; stack-specific (agent-browser tool).
- **Key content**: Core workflow (navigate, snapshot, interact, re-snapshot); command chaining; essential commands; form submission; authentication; session persistence; data extraction; parallel sessions; diffing; ref lifecycle; templates.

### orchestrating-swarms / swarm-primitives / swarm-spawning / swarm-patterns / swarm-operations
- **Type**: skill (facade + 4 specialists)
- **Path**: `.claude/skills/orchestrating-swarms/SKILL.md`, `swarm-primitives/SKILL.md`, `swarm-spawning/SKILL.md`, `swarm-patterns/SKILL.md`, `swarm-operations/SKILL.md`
- **Proposed Layer**: L2
- **Relevance**: Claude Code swarm orchestration; architecture patterns for multi-agent workflows; stack-specific (Claude Code v2.1.45).
- **Key content**: Teams, teammates, tasks, inboxes; Task vs Task+team_name+name; agent types (Explore, Plan, general-purpose, Bash, plugin agents); spawn backends (in-process, tmux, iterm2); TeamCreate, SendMessage, TeamDelete; TaskCreate, TaskUpdate, TaskList, TaskGet; patterns (parallel specialists, pipeline, swarm, research+implement, plan approval, coordinated refactoring); shutdown sequence; message formats.

### external-pattern-integrator
- **Type**: skill
- **Path**: `.claude/skills/external-pattern-integrator/SKILL.md`
- **Proposed Layer**: L2
- **Relevance**: Integrates external patterns (GSD, BMAD-METHOD) into local skills; stack research; reference examples.
- **Key content**: Phase 1 (parallel candidate mapping, tracking document); Phase 2 (contextual enhancement, workflow stage fit, external artifact recognition); Phase 3 (validation, deferred to backlog); external framework artifacts (STATE.md, ROADMAP.md, *.agent.yaml, etc.).

### universe-of-thoughts
- **Type**: skill
- **Path**: `.claude/skills/universe-of-thoughts/SKILL.md`
- **Proposed Layer**: L2 (creative reasoning)
- **Relevance**: Creative reasoning for ill-defined problems; exploratory/transformative paradigms.
- **Key content**: When to use (ambiguous goals, vast solution space, innovation); when NOT (well-defined, math, convergent); Combinational → Exploratory → Transformative; problem decomposition; evaluation criteria (feasibility, utility, novelty).

---

## ARL Meta-Layer (Observe, Identify, Probe, Accumulate, Improve)

### session-historian
- **Type**: skill
- **Path**: `.claude/skills/session-historian/SKILL.md`
- **Proposed Layer**: ARL (Observe)
- **Relevance**: Observe prior sessions; recall context when lost; search JSONL transcripts.
- **Key content**: session_query.py (list, messages, search, show, index); workflow for "I forgot what happened"; summary template; fidelity rules (read before summarizing, verbatim user messages, preserve counts, distinguish absence); index location (~/.claude/kaizen/).

### knowledge-explorer
- **Type**: skill
- **Path**: `.claude/skills/knowledge-explorer/SKILL.md`
- **Proposed Layer**: ARL (Accumulate)
- **Relevance**: Manage research/ KB; local domain knowledge with staleness (verified, next_review).
- **Key content**: list, show-template, fetch-github, add, update-append, migrate; skill-spec frontmatter (topic, category, verified, next_review); valid categories; staleness via next_review.

### refresh-research
- **Type**: skill
- **Path**: `.claude/skills/refresh-research/SKILL.md`
- **Proposed Layer**: ARL (Accumulate)
- **Relevance**: Bulk-refresh research entries; staleness detection; RT-ICA pre-flight; waves of research-curator agents.
- **Key content**: Inventory and staleness (Freshness Tracking, past review date, >6 months); scope filter (--all, --stale, --category); RT-ICA pre-flight; spawn agents in waves of 5; update README; summary report; post-actions (lint, commit).

### research-curator
- **Type**: skill
- **Path**: `.claude/skills/research-curator/SKILL.md`
- **Proposed Layer**: ARL (Accumulate)
- **Relevance**: Add/maintain research entries; single URL, batch, rerun, validate modes; stack research.
- **Key content**: Mode routing; default (single URL), batch (waves of 5), rerun (re-research), validate (script + agent fix); post-actions (README, lint, commit, push); output formats.

### research-curator/references/entry-template.md
- **Type**: reference
- **Path**: `.claude/skills/research-curator/references/entry-template.md`
- **Proposed Layer**: ARL (Accumulate)
- **Relevance**: Standard format for research entries; category selection flowchart; freshness tracking.
- **Key content**: Category selection (agent-frameworks, developer-tools, etc.); required information (identity, substance, relevance); entry file template; Freshness Tracking (Last Verified, Next Review).

### research-curator/references/validation-rules.md
- **Type**: reference
- **Path**: `.claude/skills/research-curator/references/validation-rules.md`
- **Proposed Layer**: ARL (Accumulate)
- **Relevance**: Validation for research entries; error/warning/info severity; script vs agent responsibility.
- **Key content**: section_completeness, header_fields, empty_sections (error); access_dates, freshness_tracking, statistics_currency (warning); script detects, agent fixes content.

### research-curator/references/batch-mode.md
- **Type**: reference
- **Path**: `.claude/skills/research-curator/references/batch-mode.md`
- **Proposed Layer**: ARL (Accumulate)
- **Relevance**: Wave spawning for batch research; concurrency limits.
- **Key content**: URL parsing; waves of 5; sequential waves; error handling; progress reporting.

### fact-check (ARL: Identify)
- **Type**: skill
- **Path**: `.claude/skills/fact-check/SKILL.md`
- **Proposed Layer**: ARL (Identify)
- **Relevance**: Identify deviations/hallucinations; verify claims against primary sources.
- **Key content**: (See L0 entry — also serves ARL Identify for detecting unverified/refuted claims in backlog and docs.)

### optimize-claude-md (ARL: Improve)
- **Type**: skill
- **Path**: `.claude/skills/optimize-claude-md/SKILL.md`
- **Proposed Layer**: ARL (Improve)
- **Relevance**: Improve AI-facing files; SAM + connective tissues (optimization principles, verification).
- **Key content**: (See L1 entry — also serves ARL Improve via optimization of prompts and documentation.)

### work-backlog-item close path (ARL: Improve)
- **Type**: workflow
- **Path**: `.claude/skills/work-backlog-item/SKILL.md` (Step 9)
- **Proposed Layer**: ARL (Improve)
- **Relevance**: Verification agent spawn for acceptance criteria; closes loop on completed items.
- **Key content**: Checklist verification; acceptance criteria verification agent; PASS/FAIL verdict; closing record to BACKLOG.md.

---

## Summary by Layer

| Layer | Count | Key Themes |
|-------|-------|------------|
| **L0** | 14 | SAM pipeline, RT-ICA, subagent contract, fact-check, find-cause, verify, delegate, CoVe, scientific-thinking, evidence discipline, commit conventions |
| **L1** | 8 | Agent schema, role archetypes, skill research process, quality gates, optimization principles |
| **L2** | 14 | GitHub (gh, milestones, labels, projects, issues), backlog/milestone workflow, MR changelog, daily releases, agent-browser, swarm orchestration, external pattern integration, universe-of-thoughts |
| **ARL** | 8 | session-historian (Observe), knowledge-explorer, refresh-research, research-curator (Accumulate), fact-check (Identify), optimize-claude-md, work-backlog-item close (Improve) |

---

## Cross-Cutting Concepts

- **RT-ICA**: Invoked by work-backlog-item, groom-backlog-item, refresh-research, skill-research-process (implicit), optimize-claude-md
- **CoVe**: fact-check, cove-prompt-design, optimize-claude-md
- **Evidence discipline**: fact-check, find-cause, verify
- **DONE/BLOCKED**: subagent-contract, agent-creator role archetypes, work-backlog-item (RT-ICA BLOCKED)
- **Human touchpoints**: work-backlog-item (interactive browser, AskUserQuestion), create-backlog-item, create-milestone, group-items-to-milestone, start-milestone, complete-milestone
- **Artifact conventions**: sam-definition (ARTIFACT tokens), work-backlog-item (Plan field, BACKLOG.md structure), research-curator (entry template)
- **Wave spawning**: groom-backlog-item, fact-check, research-curator, refresh-research (max 5 concurrent)
