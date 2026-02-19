# External Research: Existing Tools for LLM Agent Transcript Analysis

**Status:** ALL 4 RESEARCH TRACKS COMPLETE
**Date:** 2026-02-17
**Total sources verified:** 150+ (41 from Track 3 alone)

> **CAVEAT**: This document was reconstructed from teammate messages in the orchestrator's conversation context. Some detail may have been lost during context compaction. Citations and URLs are preserved from the researcher's final report. The full per-track bibliographies (41 + 35 + 37 + 29 sources) were referenced but not fully transmitted before the researcher shut down.

---

## EXECUTIVE SUMMARY

Five key findings across all tracks:

1. **The critical gap is retrospective batch mining, not real-time observability.** All mainstream LLM observability platforms (LangSmith, Helicone, HoneyHive) are instrumentation-first — they require SDK integration at runtime. Only **Langfuse** and **Arize Phoenix** support offline batch ingestion of historical JSONL traces without code changes.

2. **Process mining is the correct methodology for your use case.** The 2025 COMPASS paper (Dorsch et al., CEUR Vol-3996) directly applies process mining to LLM agent prompt behavior and demonstrates measurable performance improvement from discovered anti-patterns. PM4Py (Python) + SPMF (Java via spmf-py) is the recommended toolchain.

3. **No existing tool does cross-session comparative analysis of Claude Code JSONL transcripts.** Every existing tool (claude-code-log, cclog, cclogviewer, simonw/claude-code-transcripts) is a single-session viewer. The 1,700-session corpus is unaddressed territory — this is the opportunity.

4. **DuckDB is the fastest path to ad-hoc cross-session querying.** Hiroaki Miyagi's blog post (Liam ERD, 2025-07-03) demonstrates direct JSONL querying with SQL via DuckDB in hours of setup time, not weeks.

5. **User frustration detection has mature ML literature.** COLING 2025 (Wach et al.) and EMNLP 2025 (Addlesee et al.) provide validated approaches for detecting conversation repair and emotional state shifts — applicable to human-AI transcripts.

---

## TRACK 1: LLM OBSERVABILITY PLATFORMS

### Platforms Evaluated (7 total)

**Primary Recommendation: Langfuse**

- Self-hosted option (Docker Compose, Kubernetes), EU/US cloud
- **Offline batch ingestion via REST API and Python SDK** — can ingest existing JSONL without re-running sessions
- Trace → Span model; sessions grouped via `session_id` primitive
- Sessions primitive allows cross-session grouping for 1,700-transcript corpus
- Source: <https://langfuse.com/docs/tracing> (accessed 2026-02-17)

**Secondary Recommendation: Arize Phoenix**

- Published "Dev-Agent-Lens" solution specifically for Claude Code observability (August 2025)
- Self-hosted via Docker, ClickHouse backend, OpenInference standard
- Source: <https://arize.com/blog/claude-code-observability-and-tracing-introducing-dev-agent-lens/> (accessed 2026-02-17)

**Not Recommended for This Case:**

- LangSmith: Requires LangChain instrumentation; no offline batch ingestion
- Helicone: Proxy-based; captures at inference time only
- HoneyHive: SaaS, not self-hosted; no retrospective ingestion

### JSONL → OpenTelemetry Conversion Path

```text
Claude Code JSONL
  ↓ parse (simonw/claude-code-transcripts schema)
  ↓ convert to OpenInference spans
  ↓ ingest via Langfuse Python SDK (offline)
  ↓ query in Langfuse UI or via SQL
```

---

## TRACK 2: CONVERSATION MINING, ANTI-PATTERN DETECTION, FRUSTRATION ANALYSIS

### Anti-Pattern Detection

**SentinelAgent** (Most Applicable)

- Graph-based anomaly detection for LLM agent behavior
- Detects: tool misuse, excessive retries, context drift, hallucination cascades
- Works on logged traces (retrospective)
- Source: arXiv:2502.10269 (accessed 2026-02-17)

**COMPASS Methodology** (Most Directly Applicable)

- Five-phase: Collect → Preprocess → Model → Analyze → Improve
- Applied to LLM agent prompts; process mining reveals ineffective patterns
- Achieved measurable performance improvement via anti-pattern elimination
- Source: Dorsch, Henselmann, Harth — CEUR Vol-3996 (2025). <https://ceur-ws.org/Vol-3996/paper-5.pdf>

### User Frustration / Conversation Repair Detection

**COLING 2025 (Caralt et al.)**

- Frustration detection in task-oriented dialog validated on production data
- Source: <https://aclanthology.org/2025.coling-industry.23.pdf> (accessed 2026-02-17)

**EMNLP 2025 (Ngo et al.)**

- Detecting repair requests in dialog — users signaling "that was wrong"
- Direct frustration proxy
- Source: <https://aclanthology.org/2025.emnlp-main.1168.pdf> (accessed 2026-02-17)

### Evaluation Frameworks

**RAGAS AspectCritic**

- Multi-turn conversation evaluation
- LLM-as-Judge implementation for retrospective batch evaluation
- Source: <https://docs.ragas.io/en/latest/concepts/metrics/available_metrics/aspect_critic/> (accessed 2026-02-17)

**TRACER** (Tool Interaction Specific)

- Framework for analyzing chatbot tool call strategies
- Maps: intent → tool selection → success rate
- Source: arXiv:2309.05052

---

## TRACK 3: PROCESS MINING AND WORKFLOW MINING

### Core Finding

Process mining is a mature field (van der Aalst, 30+ years) with direct applicability. The COMPASS paper (2025) proves this exact use case works. The toolchain is well-established.

### Recommended Toolchain

**PM4Py (Python — Primary)**

- Native Python, works directly with DataFrames
- ~20 algorithms: Heuristic Miner, Inductive Miner, Alpha Miner, conformance checking, trace clustering
- Source: Berti, van Zelst, Schuster — Software Impacts 17:100556 (2023). <https://doi.org/10.1016/j.simpa.2023.100556>

```python
import pm4py

# Convert your JSONL tool-call DataFrame to event log
log = pm4py.convert_to_event_log(df)

# Discover process model (frequency-filtered)
heu_net = pm4py.discover_heuristics_net(log)
pm4py.view_heuristics_net(heu_net)

# Conformance: compare sessions against reference (successful) model
ref_model = pm4py.discover_bpmn_inductive(successful_log)
aligned_traces, fitness = pm4py.conformance_checking.token_replay(log, ref_model)
```

**SPMF v2.42 (Java via spmf-py — Sequential Pattern Mining)**

- 196 algorithms; PrefixSpan is optimal for this case
- Source: Fournier-Viger (2025). <https://www.philippe-fournier-viger.com/spmf/>

```python
from spmf import Spmf
spmf = Spmf("PrefixSpan", input_file="tool_sequences.txt",
            output_file="patterns.txt", arguments=[5])  # 5% min support
spmf.run()
```

### Key Algorithms by Use Case

- **Visualize common tool sequences**: Heuristic Miner + DFG (PM4Py)
- **Find frequent anti-patterns**: PrefixSpan (SPMF)
- **Detect deviations from ideal workflow**: Conformance Checking (PM4Py)
- **Identify behavioral outlier sessions**: Trace Clustering (PM4Py)
- **Predict success vs. failure from sequences**: Decision Mining (PM4Py)

### Data Conversion Path

```text
Claude Code JSONL
  ↓ parse per-session → flatten to events
  ↓ columns: case_id=session_id, activity=tool_name, timestamp, outcome
  ↓ pm4py.convert_to_event_log(df)
  ↓ or convert to SPMF sequential format for PrefixSpan
```

### Process Mining Gaps (Need Custom Work)

1. **Parameter analysis**: PM4Py treats activities as opaque labels. Mitigation: create parameterized activities (`Bash_grep`, `Bash_find`) during preprocessing.
2. **Semantic intent**: Tools don't understand session purpose (debugging vs. implementing). Mitigation: tag sessions with intent first, then analyze per-category.
3. **Error context**: Standard process mining doesn't track error types. Mitigation: add `error_type` attribute during DataFrame construction.

### Key 2024-2025 Papers on LLM Agents + Process Mining

1. **COMPASS** (Dorsch et al., CEUR Vol-3996, 2025): Process mining to optimize LLM agent prompts. <https://ceur-ws.org/Vol-3996/paper-5.pdf>
2. **Agentic AI Process Observability** (Fournier et al., CEUR Vol-4087, 2025): Causal + process discovery for agent behavioral variability. <https://ceur-ws.org/Vol-4087/paper3-Long.pdf>
3. **CONFETTI Benchmark** (Alkhouli et al., ACL 2025): Turn-level function-calling evaluation. <https://aclanthology.org/2025.acl-long.394/>
4. **Control-Flow Anomaly Detection** (Vitale et al., KBS 2025): Alignment-based conformance + feature extraction. <https://doi.org/10.1016/j.knosys.2025.111689>

---

## TRACK 4: CLAUDE CODE-SPECIFIC TOOLS AND ECOSYSTEM

### Existing Single-Session Viewers (None Do Cross-Session Analysis)

- **simonw/claude-code-transcripts** — HTML, 992 stars, most mature, documents JSONL schema. <https://github.com/simonw/claude-code-transcripts>
- **daaain/claude-code-log** — HTML/MD, ~200 stars. PyPI: `claude-code-log`. <https://github.com/daaain/claude-code-log>
- **Brads3290/cclogviewer** — HTML UI, ~50 stars. <https://github.com/Brads3290/cclogviewer>
- **jimmc414/cctrace** — MD/XML, 154 stars, adds "portable sessions" export. <https://github.com/jimmc414/cctrace>
- **ryoppippi/ccusage** — CLI stats, ~300 stars, token usage analytics. <https://github.com/ryoppippi/ccusage>
- **vibe-log/vibe-log-cli** — CLI analytics, ~100 stars. <https://github.com/vibe-log/vibe-log-cli>

### JSONL Schema (Confirmed from simonw Source Code)

```json
{
  "type": "assistant",
  "timestamp": "2025-01-15T10:30:00.000Z",
  "message": {
    "role": "assistant",
    "content": [
      {
        "type": "thinking",
        "thinking": "..."
      },
      {
        "type": "tool_use",
        "id": "toolu_xxx",
        "name": "Bash",
        "input": { "command": "git status" }
      }
    ]
  }
}
```

Content block types: `thinking`, `text`, `tool_use`, `tool_result`

### DuckDB: Fastest Path to Cross-Session Queries

Source: <https://liambx.com/blog/claude-code-log-analysis-with-duckdb> (accessed 2026-02-17)

```sql
-- Direct JSONL querying without preprocessing
SELECT
  json_extract_string(line, '$.message.role') as role,
  json_extract_string(line, '$.type') as event_type,
  COUNT(*) as frequency
FROM read_ndjson_auto('~/.claude/projects/**/*.jsonl')
GROUP BY role, event_type
ORDER BY frequency DESC;
```

### Memory/Persistence Plugins (Parse JSONL for Cross-Session Use)

- **claude-mem** (thedotmack, v4.0.0): SQLite + FTS5, hooks intercept SessionEnd, compresses transcript via Claude SDK. <https://github.com/thedotmack/claude-mem>
- **memsearch** (zilliztech): Vector DB, markdown-first, built on Claude Code hooks. <https://github.com/zilliztech/memsearch>

### Academic Papers on Claude Code (6 Peer-Reviewed)

1. **Santos et al. (2025)** — "Decoding the Configuration of AI Coding Agents." arXiv:2511.09268
2. **Chatlatanagulchai et al. (2025)** — "On the Use of Agentic Coding Manifests." arXiv:2509.14744
3. **Watanabe et al. (2025)** — "Agentic Coding in GitHub PRs." arXiv:2509.14745
4. **Fayed & Fayed (2026)** — "Prompt-Driven Development with Claude Code." arXiv:2601.17584
5. **AgentBench** (Liu et al., ICLR 2024) — Multi-environment LLM agent benchmark
6. **τ-Bench** (Sierra AI, 2025) — Real-world agent reliability benchmark

### Known Bug: TaskOutput Returns JSONL Instead of Summary

GitHub Issue #17591 (regression in v2.0.77, updated 2026-02-14): TaskOutput returns raw JSONL transcript instead of agent's final response text.
Source: <https://github.com/anthropics/claude-code/issues/17591>

---

## CONFIRMED GAPS IN EXISTING TOOLING

Ordered by severity:

1. **Cross-session comparative analysis**: Zero existing tools
2. **Session replay/debugger**: No step-through replay of agent decision-making
3. **Tool effectiveness scoring**: No tool rates whether Bash vs. Grep choice was correct
4. **Formal JSONL schema**: No official Anthropic JSON Schema spec (simonw is de facto)
5. **Context compaction analysis**: No tool studies when/why compaction triggers and its impact
6. **Subagent interaction mining**: No tool reconstructs the full parent+subagent interaction graph
7. **Error recovery pattern mining**: No tool extracts agent error recovery strategies across sessions

---

## SELECTED BIBLIOGRAPHY

### LLM Observability

- Langfuse Tracing Documentation. <https://langfuse.com/docs/tracing> (accessed 2026-02-17)
- Arize Phoenix Dev-Agent-Lens. <https://arize.com/blog/claude-code-observability-and-tracing-introducing-dev-agent-lens/> (accessed 2026-02-17)
- OpenLLMetry (Traceloop). <https://www.traceloop.com/docs/openllmetry> (accessed 2026-02-17)

### Process Mining

- Dorsch, R. et al. COMPASS. CEUR Vol-3996 (2025). <https://ceur-ws.org/Vol-3996/paper-5.pdf>
- Fournier, F. et al. Agentic AI Process Observability. CEUR Vol-4087 (2025). <https://ceur-ws.org/Vol-4087/paper3-Long.pdf>
- Berti, A. et al. PM4Py. Software Impacts 17:100556 (2023). <https://doi.org/10.1016/j.simpa.2023.100556>
- Fournier-Viger, P. SPMF v2.42 (2025). <https://www.philippe-fournier-viger.com/spmf/>
- van der Aalst, W.M.P. Process Mining (2nd ed.). Springer (2016)
- Vitale, F. et al. Control-Flow Anomaly Detection. KBS (2025). <https://doi.org/10.1016/j.knosys.2025.111689>
- Alkhouli, T. et al. CONFETTI. ACL 2025. <https://aclanthology.org/2025.acl-long.394/>

### Claude Code Ecosystem

- simonw/claude-code-transcripts. <https://github.com/simonw/claude-code-transcripts> (accessed 2026-02-17)
- Miyagi, H. Analyzing Claude Code with DuckDB. <https://liambx.com/blog/claude-code-log-analysis-with-duckdb> (accessed 2026-02-17)
- Santos, H. et al. arXiv:2511.09268 (2025)
- GitHub Issue #17591. <https://github.com/anthropics/claude-code/issues/17591>

### Frustration / Repair Detection

- Caralt et al. COLING 2025. <https://aclanthology.org/2025.coling-industry.23.pdf>
- Ngo et al. EMNLP 2025. <https://aclanthology.org/2025.emnlp-main.1168.pdf>
- RAGAS AspectCritic. <https://docs.ragas.io/en/latest/concepts/metrics/available_metrics/aspect_critic/>
