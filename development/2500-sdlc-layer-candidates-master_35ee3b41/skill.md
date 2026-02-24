# SDLC Layer Separation Architecture — Master Candidate List

**Generated**: 2026-02-23  
**Source**: Full read of plugins/, .claude/skills/, .claude/agents/ by 3 parallel sub-agents  
**Purpose**: Comprehensive review list for user to add specific items to the plan

---

## How to Use This List

1. **Review each section** by layer (L0, L1, L2, ARL)
2. **Mark candidates** for inclusion: add to plan, defer, or note as "already covered"
3. **Extract** specific concepts, diagrams, workflows, references, KB articles, paths

**Detail reports** (full per-item formats):
- `sdlc-layer-candidates-from-plugins.md` — plugins
- `sdlc-layer-candidates.md` — skills
- `sdlc-layer-agent-candidates.md` — agents

---

## Layer 0 (SDLC-Agnostic) — Candidates for Review

### From Plugins

| Item | Path | Key Content |
|------|------|-------------|
| development-harness | `plugins/development-harness/CLAUDE.md` | SAM 7-stage pipeline, artifact conventions, human touchpoints, Voltron composition |
| default-development-flow | `plugins/development-harness/.../default-development-flow.md` | Mermaid flowchart, Gate 1/2, loop limits (3 NEEDS_WORK, 2 NOT_CERTIFIED) |
| human-touchpoint-model | `plugins/development-harness/.../human-touchpoint-model.md` | Bound/Unbound/Mixed constraints, risk assessment, escalation triggers |
| artifact-conventions | `plugins/development-harness/.../artifact-conventions.md` | ARTIFACT:{TYPE}({ID}), required sections, `.planning/harness/` |
| role-resolution-protocol | `plugins/development-harness/.../role-resolution-protocol.md` | Abstract roles → agents, project detection, manifest search |
| language-manifest-schema | `plugins/development-harness/.../language-manifest-schema.md` | Role Fulfillment, Quality Gates, Project Detection, Flow Override |
| Workflow stages (7) | `plugins/development-harness/skills/workflows/*/SKILL.md` | Discovery, Planning, Context Integration, Task Decomposition, Execution, Forensic Review, Final Verification |
| planner-rt-ica | `plugins/development-harness/skills/planner-rt-ica/` | APPROVED-FOR-PLANNING / APPROVED-WITH-GAPS / BLOCKED-FOR-PLANNING |
| validation-protocol, clear-cove-task-design, generate-task | `plugins/development-harness/skills/` | Validation checklists, CoVe task design, task file format |
| verification-gate | `plugins/verification-gate/skills/verification-gate/SKILL.md` | 4-checkpoint pre-action verification, hypothesis-action alignment |
| orchestrator-discipline | `plugins/orchestrator-discipline/` | Read constraints, Investigation Escalation Anti-Pattern, Agent Output Polling Anti-Pattern |
| agent-orchestration | `plugins/agent-orchestration/skills/agent-orchestration/SKILL.md` | WHERE-WHAT-WHY, pre-gathering anti-pattern, delegation template |
| how-to-delegate | `plugins/agent-orchestration/skills/how-to-delegate/SKILL.md` | 10-step delegation worksheet, pre-flight verification |
| hallucination-detector | `plugins/hallucination-detector/` | 5 triggers (Speculation, Causality, Pseudo-Quantification, Completeness, Micromanagement) |
| the-rewrite-room/output-contracts | `plugins/the-rewrite-room/.../output-contracts.md` | STATUS block (DONE/BLOCKED/FAILED), artifact listing |
| the-rewrite-room/fidelity-rules | `plugins/the-rewrite-room/.../fidelity-rules.md` | Read before abstracting, preserve counts, distinguish absence |
| plugin-creator/arl | `plugins/plugin-creator/skills/arl/SKILL.md` | ARL theory, 10 gates (R1-R10), HOOTL, scope-feasibility matrix |

### From Skills

| Item | Path | Key Content |
|------|------|-------------|
| work-backlog-item | `.claude/skills/work-backlog-item/SKILL.md` | BACKLOG → SAM bridge, --auto, close/resolve paths |
| sam-definition | `.claude/skills/work-backlog-item/references/sam-definition.md` | Canonical SAM, artifact flow, escalation rules |
| github-integration | `.claude/skills/work-backlog-item/references/github-integration.md` | BACKLOG ↔ Issue mapping, setup-github |
| validation-plan | `.claude/skills/work-backlog-item/references/validation-plan.md` | V1–V6 verification commands |
| groom-backlog-item | `.claude/skills/groom-backlog-item/SKILL.md` | fact-check → RT-ICA → groomer flow |
| rt-ica | `.claude/skills/rt-ica/SKILL.md` | Goal reconstruction, AVAILABLE/DERIVABLE/MISSING, BLOCK |
| subagent-contract | `.claude/skills/subagent-contract/SKILL.md` | DONE/BLOCKED, scope discipline, supervisor protocol |
| fact-check | `.claude/skills/fact-check/SKILL.md` | Primary-source verification, wave spawning, VERIFIED/REFUTED/INCONCLUSIVE |
| find-cause | `.claude/skills/find-cause/SKILL.md` | Evidence chain, reproduce, observe, source read |
| verify | `.claude/skills/verify/SKILL.md` | Completion checklist, WORKS/FIXED checks |
| delegate | `.claude/skills/delegate/SKILL.md` | OBSERVATIONS, DEFINITION OF SUCCESS, CONTEXT, YOUR TASK |
| cove-prompt-design | `.claude/skills/cove-prompt-design/SKILL.md` | Chain of Verification, verification questions |
| scientific-thinking | `.claude/skills/scientific-thinking/SKILL.md` | Hypothesis → Prediction → Experiment → Conclusion |
| commit-staged | `.claude/skills/commit-staged/SKILL.md` | Conventional commits, feat/fix/docs/scope |

### From Agents

| Item | Path | Key Content |
|------|------|-------------|
| backlog-item-groomer | `.claude/agents/backlog-item-groomer.md` | Context manifest, RT-ICA procedure, dependency graph |
| fact-checker | `.claude/agents/fact-checker.md` | CoVe, structured verdict, primary sources only |
| topic-specialist | `.claude/agents/topic-specialist.md` | Evidence discipline, skill update/create with verified findings |
| context-gathering | `.claude/agents/context-gathering.md` | Context manifest, "How It Currently Works" narrative |
| context-refinement | `.claude/agents/context-refinement.md` | Drift detection, "Discovered During Implementation" |
| logging | `.claude/agents/logging.md` | Work Log consolidation, chronological order |
| code-review | `.claude/agents/code-review.md` | LLM-slop detection, Critical/Warning/Suggestion |

---

## Layer 1 (Language-Specific) — Candidates for Review

### From Plugins

| Item | Path | Key Content |
|------|------|-------------|
| language-manifest-template | `plugins/development-harness/templates/language-manifest-template.md` | Skeleton for Python, TypeScript, Rust, Go |
| python3-development | `plugins/python3-development/skills/python3-development/SKILL.md` | ROLE_TYPE, Linting Discovery, Format-First, quality gates |
| python-development-orchestration | `plugins/python3-development/.../python-development-orchestration.md` | TDD, Feature Addition, Code Review workflows |
| perl-development | `plugins/perl-development/` | perl-script-developer, perl-script-auditor, perl-cli-architect |
| bash-development | `plugins/bash-development/` | bash-script-developer, bash-script-auditor, version-specific features |

### From Skills

| Item | Path | Key Content |
|------|------|-------------|
| agent-creator | `.claude/skills/agent-creator/SKILL.md` | Abstract roles, Phase 1–7 workflow |
| agent-schema | `.claude/skills/agent-creator/references/agent-schema.md` | Frontmatter spec, language manifest for agents |
| agent-templates | `.claude/skills/agent-creator/references/agent-templates.md` | Role archetypes: Researcher, Planner, Coder, Tester, Reviewer |
| agent-examples | `.claude/skills/agent-creator/references/agent-examples.md` | Real-world agent implementations |
| skill-research-process | `.claude/skills/skill-research-process/SKILL.md` | Quality gates, anti-hallucination checkpoint |
| optimize-claude-md | `.claude/skills/optimize-claude-md/SKILL.md` | 8 optimization principles, RT-ICA/CoVe integration |

### From Agents

| Item | Path | Key Content |
|------|------|-------------|
| javascript-pro | `.claude/agents/javascript-pro.md` | ES2023+, Node.js, quality checklist, project detection |
| typescript-pro | `.claude/agents/typescript-pro.md` | TypeScript 5.0+, strict mode, type coverage |
| c-systems-programmer | `.claude/agents/c-systems-programmer.md` | C99/C11, POSIX, Valgrind, GDB |

---

## Layer 2 (Stack/Goal-Specific) — Candidates for Review

### From Plugins

| Item | Path | Key Content |
|------|------|-------------|
| the-rewrite-room | `plugins/the-rewrite-room/skills/the-rewrite-room/SKILL.md` | Workflow router, 7 workflows, adapters |
| the-rewrite-room/workflows | `plugins/the-rewrite-room/.../workflows/*.md` | drift-audit, documentation-sync, authoring, prompt-optimization, etc. |
| summarizer | `plugins/summarizer/skills/summarizer/SKILL.md` | Format routing, fidelity rules, team coordination |
| holistic-linting | `plugins/holistic-linting/` | Lint orchestration, root-cause resolver |
| gitlab-skill | `plugins/gitlab-skill/` | GitLab CI, GLFM, gitlab-ci-local |
| fastmcp-creator | `plugins/fastmcp-creator/` | MCP server creation, FastMCP 3.x |
| dasel | `plugins/dasel/` | Structured data query/transform |
| clang-format | `plugins/clang-format/` | C/C++ formatting config |
| xdg-base-directory | `plugins/xdg-base-directory/` | Cross-platform config storage |
| plugin-creator | `plugins/plugin-creator/` | Plugin creation, agent-creator, refactor workflows |
| development-harness agents | `plugins/development-harness/agents/*.md` | swarm-task-planner, plan-validator, feature-researcher, codebase-analyzer, etc. |

### From Skills

| Item | Path | Key Content |
|------|------|-------------|
| gh | `.claude/skills/gh/SKILL.md` | GitHub CLI, -R flag, proxy environments |
| gh/milestones, labels, projects-v2, issue-stories | `.claude/skills/gh/references/` | Milestone CRUD, label taxonomy, Projects V2, issue lifecycle |
| create-backlog-item | `.claude/skills/create-backlog-item/SKILL.md` | Guided/quick/--auto modes |
| create-milestone, group-items-to-milestone, start-milestone, complete-milestone | `.claude/skills/` | Milestone workflow |
| create-merge-request-changelog | `.claude/skills/create-merge-request-changelog/SKILL.md` | Domain-based change categorization |
| daily-releases | `.claude/skills/daily-releases/SKILL.md` | AI changelogs per day |
| agent-browser | `.claude/skills/agent-browser/SKILL.md` | Browser automation CLI |
| orchestrating-swarms / swarm-primitives / swarm-spawning / swarm-patterns / swarm-operations | `.claude/skills/` | Teams, tasks, inboxes, patterns |
| external-pattern-integrator | `.claude/skills/external-pattern-integrator/SKILL.md` | GSD, BMAD-METHOD integration |
| universe-of-thoughts | `.claude/skills/universe-of-thoughts/SKILL.md` | Creative reasoning for ill-defined problems |

### From Agents

| Item | Path | Key Content |
|------|------|-------------|
| plugin-assessor | `.claude/agents/plugin-assessor.md` | Plugin schema validation, orphan classification |
| plugin-docs-writer | `.claude/agents/plugin-docs-writer.md` | User-facing README generation |
| research-curator | `.claude/agents/research-curator.md` | Research entry creation, categories |
| research-context-agent | `.claude/agents/research-context-agent.md` | Integration opportunities, cross-references |

---

## ARL Meta-Layer — Candidates for Review

### From Plugins

| Item | Path | Key Content |
|------|------|-------------|
| agentskill-kaizen | `plugins/agentskill-kaizen/` | Observe: transcript analysis; Improve: kaizen-improvement |
| transcript-analysis | `plugins/agentskill-kaizen/skills/transcript-analysis/SKILL.md` | 10 dimensions, DuckDB, process mining |
| kaizen-improvement | `plugins/agentskill-kaizen/skills/kaizen-improvement/SKILL.md` | Hook generation, agent/skill patches |
| hallucination-detector | `plugins/hallucination-detector/` | Identify: speculation-as-diagnosis, invented causality |
| plugin-creator/arl | `plugins/plugin-creator/skills/arl/SKILL.md` | ARL theory, Layer 3 Observation |

### From Skills

| Item | Path | Key Content |
|------|------|-------------|
| session-historian | `.claude/skills/session-historian/SKILL.md` | Observe: session recall, transcript search |
| knowledge-explorer | `.claude/skills/knowledge-explorer/SKILL.md` | Accumulate: research/ KB, staleness |
| refresh-research | `.claude/skills/refresh-research/SKILL.md` | Accumulate: bulk refresh, waves of 5 |
| research-curator | `.claude/skills/research-curator/SKILL.md` | Accumulate: add/maintain research entries |
| research-curator/entry-template, validation-rules, batch-mode | `.claude/skills/research-curator/references/` | Entry format, freshness tracking |
| fact-check (ARL: Identify) | `.claude/skills/fact-check/SKILL.md` | Deviations/hallucinations detection |
| optimize-claude-md (ARL: Improve) | `.claude/skills/optimize-claude-md/SKILL.md` | Improve AI-facing files |
| work-backlog-item close path | `.claude/skills/work-backlog-item/SKILL.md` Step 9 | Verification agent, acceptance criteria |

### From Agents

| Item | Path | Key Content |
|------|------|-------------|
| doc-drift-auditor | `.claude/agents/doc-drift-auditor.md` | Identify: doc vs implementation drift |
| code-review (ARL overlap) | `.claude/agents/code-review.md` | Identify: LLM-slop, hallucination patterns |
| context-refinement (ARL overlap) | `.claude/agents/context-refinement.md` | Accumulate: institutional knowledge |
| logging (ARL overlap) | `.claude/agents/logging.md` | Observe: work session consolidation |
| topic-specialist (ARL overlap) | `.claude/agents/topic-specialist.md` | Improve: skill update/create with verified findings |
| research-context-agent (ARL overlap) | `.claude/agents/research-context-agent.md` | Accumulate: integration opportunities |

---

## Cross-Cutting Concepts (Add to Plan)

| Concept | Sources |
|---------|---------|
| **RT-ICA** | work-backlog-item, groom-backlog-item, refresh-research, optimize-claude-md, backlog-item-groomer |
| **CoVe** | fact-check, cove-prompt-design, optimize-claude-md, fact-checker, topic-specialist |
| **Evidence discipline** | fact-check, find-cause, verify, fact-checker, topic-specialist |
| **DONE/BLOCKED** | subagent-contract, agent-creator role archetypes, work-backlog-item |
| **Human touchpoints** | work-backlog-item, create-backlog-item, create-milestone, group-items-to-milestone, start-milestone, complete-milestone |
| **Artifact conventions** | sam-definition, work-backlog-item, artifact-conventions, output-contracts |
| **Wave spawning** | groom-backlog-item, fact-check, research-curator, refresh-research (max 5 concurrent) |

---

## Plugins Not Yet Mapped (Minimal SDLC Relevance)

- brainstorming-skill, commitlint, conventional-commits, litellm, llamafile, prompt-optimization-claude-45, uv

---

## Next Steps

1. **Add to plan**: Copy specific items from this list into the SDLC Layer Separation Architecture plan (Layer 0 hub, Layer 1 formalization, Layer 2 schema, ARL meta-layer).
2. **Extract mermaid diagrams**: default-development-flow, role-resolution-protocol, others.
3. **Extract methodologies**: orchestrator-discipline, agent-orchestration, verification-gate.
4. **Map research/**: Run fourth sub-agent to inventory research/ and .claude/docs/ for KB articles and reference docs.
