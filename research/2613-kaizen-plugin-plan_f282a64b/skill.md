# Plugin Plan: agentskill-kaizen

**Status:** Phase 3 — Detailed Design (in progress)
**Date:** 2026-02-18

## Phase 1: Discovery

**Purpose:** Analyze Claude Code agent session transcripts to identify inefficiencies, anti-patterns, repeated mistakes, missing tooling opportunities, and user frustration signals — enabling continuous improvement (kaizen) of agents, skills, and workflows.

**Problem:** Agents make the same mistakes across sessions, repeat multi-step workflows that could be scripts, and frustrate users by ignoring instructions. No systematic way to detect these patterns. Data exists in ~1,700+ local JSONL transcripts. No existing tool does cross-session analysis.

**Target users:** Plugin developers and power users maintaining agent/skill definitions.

**Plugin type:** Analysis toolkit

### Analysis Dimensions (10)

1. Tool misuse / inefficiency (Bash ls instead of Glob, etc.)
2. Repeated errors across sessions
3. User frustration signals (corrections, interrupts, "no don't do that")
4. Missing tooling opportunities (repeated manual workflows → scripts/skills)
5. Subagent delegation patterns (410 Task calls, 30+ agent types)
6. Shortest path to outcome — where are steps being wasted?
7. Red herrings — distractions chased across sessions that never lead to the solution
8. System process interruptions — hooks, permission prompts, compaction interrupting correct paths
9. Missing hooks — behavioral corrections done manually that should be automated
10. MCP-based DuckDB querying — MotherDuck MCP for direct SQL access to transcripts

### Research Documents

- `.claude/kaizen-data-analysis.md` — JSONL schema, signal catalog, pipeline architecture (764 lines)
- `.claude/kaizen-external-research.md` — existing tools, frameworks, academic papers, confirmed gaps
- `.claude/kaizen-references.md` — hooks API, headless mode, MotherDuck MCP, expanded dimensions

### Key Research Findings

- ~500MB corpus, ~57k JSONL records, 13 record types, 687 subagent transcripts
- 593 Bash file-op violations already quantified (215× ls, 182× grep, 105× head/tail/cat)
- 106 user interrupts, ~43 real correction utterances
- PM4Py (process mining) + SPMF (sequential pattern mining) are the proven toolchain
- COMPASS paper (2025) validated process mining on LLM agent traces
- No existing tool does cross-session analysis — confirmed gap

---

## Phase 2: Component Plan

| Component Type | Count | Purpose |
|----------------|-------|---------|
| Skills | 2 | `transcript-analysis` (methodology, signal catalog), `kaizen-improvement` (turning findings into outputs) |
| Commands | 4 | `/kaizen:analyze`, `/kaizen:explore`, `/kaizen:report`, `/kaizen:generate-hooks` |
| Agents | 2 | `transcript-analyst`, `improvement-generator` |
| Hooks | 0 | Deferred — start with batch analysis, add real-time capture later |
| MCP | 2 | MotherDuck DuckDB (SQL querying), custom FastMCP (process mining, pattern mining, frustration detection) |
| Scripts | 3-4 | Supporting the custom MCP server |

### Architecture Flow

```text
/kaizen:analyze (autonomous pipeline)
  → transcript-analyst agent spawns
  → queries JSONL via DuckDB MCP (execute_query)
  → calls custom MCP tools for deeper analysis (process mining, pattern mining)
  → writes findings to .planning/kaizen/

/kaizen:explore (interactive mode)
  → transcript-analyst agent in interactive mode
  → user steers investigation
  → same tools available

/kaizen:report
  → reads existing analysis in .planning/kaizen/
  → generates summary report

/kaizen:generate-hooks
  → reads analysis findings
  → improvement-generator agent produces hook proposals
  → default: draft proposals with rationale
  → --install flag: writes directly to settings
```

---

## Phase 3: Detailed Design Decisions

### Resolved via User Input

1. **Analysis scope**: User chooses at runtime via `--project` flag, defaults to current project
2. **Output format**: Markdown files in `.planning/kaizen/`
3. **Hook generation**: Draft proposals by default, `--install` flag writes to settings
4. **Custom MCP language**: Python (FastMCP) — natural fit for PM4Py, SPMF, pandas
5. **Agent mode**: `/kaizen:analyze` runs autonomous pipeline, `/kaizen:explore` enables interactive mode
6. **SessionEnd hook**: Deferred — start with batch analysis, add real-time metrics later
7. **Query approach**: Dynamic query generation — agent writes SQL on the fly
8. **Patch format**: Instruction sets for other agents — follows delegation pattern (outcome-focused prompts for contextual-ai-documentation-optimizer, skill-creator, etc.)

### Component Specifications

#### Skill: transcript-analysis

- **Triggers**: "analyze transcripts", "session analysis", "find anti-patterns", "kaizen"
- **Content**: JSONL schema reference, signal catalog with field paths, extraction methodology
- **References**: SQL query patterns, PM4Py usage, frustration signal taxonomy
- **No utility scripts** — analysis runs through MCP tools

#### Skill: kaizen-improvement

- **Triggers**: "generate hooks from findings", "improve agent", "fix anti-pattern"
- **Content**: How to translate findings into actionable outputs (hooks, agent patches, skill improvements)
- **References**: Hook API patterns, agent frontmatter format, delegation prompt templates

#### Command: /kaizen:analyze

- **Arguments**: `--project <name>` (scope), `--dimensions <list>` (which analysis dimensions)
- **Tools**: DuckDB MCP, custom kaizen MCP, Read, Write, Glob
- **Output**: Writes analysis report to `.planning/kaizen/analysis-{date}.md`

#### Command: /kaizen:explore

- **Arguments**: `--project <name>` (scope)
- **Tools**: Same as analyze
- **Mode**: Interactive — presents initial findings, user steers deeper

#### Command: /kaizen:report

- **Arguments**: `--latest` (most recent analysis), `--all` (aggregate)
- **Tools**: Read, Write
- **Output**: Summary report in `.planning/kaizen/report-{date}.md`

#### Command: /kaizen:generate-hooks

- **Arguments**: `--install` (write to settings), `--from <analysis-file>` (source findings)
- **Tools**: Read, Write, Edit
- **Output**: Hook proposals in `.planning/kaizen/hooks/` or directly to settings with --install

#### Agent: transcript-analyst

- **Model**: sonnet
- **Tools**: DuckDB MCP (execute_query), custom kaizen MCP, Read, Glob, Grep, Bash
- **Skills**: transcript-analysis
- **Trigger**: Spawned by /kaizen:analyze and /kaizen:explore commands
- **Output**: Structured analysis findings written to .planning/kaizen/

#### Agent: improvement-generator

- **Model**: sonnet
- **Tools**: Read, Write, Glob, Grep
- **Skills**: kaizen-improvement
- **Trigger**: Spawned by /kaizen:generate-hooks or after analysis completes
- **Output**: Outcome-focused instruction sets delegatable to other agents (optimizer, skill-creator)

#### MCP: MotherDuck DuckDB

- **Config**: `uvx mcp-server-motherduck --db-path :memory: --read-write`
- **Purpose**: Direct SQL querying of JSONL transcript files
- **Tools exposed**: execute_query, list_tables, list_columns

#### MCP: Custom Kaizen Analysis Server (FastMCP)

- **Language**: Python with FastMCP
- **Tools to expose**:
  - `discover_process_model` — PM4Py Heuristic Miner on tool-call sequences
  - `check_conformance` — PM4Py conformance checking against reference model
  - `find_frequent_patterns` — SPMF PrefixSpan on tool-call sequences
  - `detect_frustration_signals` — NLP extraction from user turns
  - `cluster_sessions` — PM4Py trace clustering by behavioral similarity
  - `extract_tool_sequences` — JSONL → ordered tool-call arrays per session
- **Dependencies**: pm4py, spmf-py, pandas, fastmcp

---

## Phase 4-8: Not yet started

- Phase 4: Plugin structure creation
- Phase 5: Component implementation
- Phase 6: Validation & quality check
- Phase 7: Testing & verification
- Phase 8: Documentation & next steps
