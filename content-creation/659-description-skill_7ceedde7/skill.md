---
description: Route summarization requests to the correct methodology and enforce fidelity rules. Activates on summarize, tl;dr, give me the highlights, what's important in this, break down this, what does this code do, explain this file, describe this image, read and summarize. Routes files, URLs, images, and multi-source content to type-specific strategies. Enforces anti-hallucination rules — read before summarizing, extract before abstracting, preserve counts, distinguish absence from nonexistence, prevent lossy re-summarization chains.
---
# Summarizer

Route to the correct summarization methodology and enforce fidelity rules across all summarization operations.

## Format Selection

When the user specifies an output format, load the corresponding template before delegating to a domain skill or agent. If no format is specified, default to `structured`.

| User Signal | format_id | Template |
|-------------|-----------|----------|
| (no format specified) | `structured` | [structured.md](./templates/structured.md) |
| "bullet points", "key points", "quick bullets" | `bullets` | [bullets.md](./templates/bullets.md) |
| "tl;dr", "one-liner", "in a nutshell" | `tldr` | [tldr.md](./templates/tldr.md) |
| "json", "machine-readable", "structured data" | `json` | [json.md](./templates/json.md) |
| "table", "tabular", "grid format" | `table` | [table.md](./templates/table.md) |
| "outline", "table of contents", "hierarchical" | `outline` | [outline.md](./templates/outline.md) |

**Loading instruction**: Read `$SKILL_DIR/templates/{format_id}.md` to obtain the schema, example, and fidelity constraints for the selected format. Pass the `format` parameter to the delegated agent or apply the template directly when summarizing inline.

## Decision Tree

When the model needs to summarize content, follow this decision tree to select the correct approach:

```text
INPUT RECEIVED
  │
  ├─ Is it a FILE path?
  │   ├─ Yes → Load file-summarization skill
  │   │         Run file_metrics.py to assess size
  │   │         Select strategy based on size thresholds
  │   └─ No ↓
  │
  ├─ Is it a URL?
  │   ├─ Yes → Load url-summarization skill
  │   │         Fetch content first, then summarize
  │   └─ No ↓
  │
  ├─ Is it an IMAGE (path to .png, .jpg, .gif, .svg, .webp, screenshot)?
  │   ├─ Yes → Load image-summarization skill
  │   │         Read image with Read tool (multimodal)
  │   └─ No ↓
  │
  ├─ Is it MULTIPLE sources (list of files, URLs, or mixed)?
  │   ├─ Yes → Count sources
  │   │   ├─ 3+ sources AND Teammate tool available?
  │   │   │   ├─ Yes → Spawn summarizer team (see Team Coordination below)
  │   │   │   │         Each teammate summarizes one source
  │   │   │   │         Teammates cross-check findings via messaging
  │   │   │   │         Leader synthesizes with multi-source-synthesis skill
  │   │   │   └─ No ↓
  │   │   └─ Summarize each source individually (subagents or sequential)
  │   │             Then load multi-source-synthesis skill
  │   │             Combine with deduplication and attribution
  │   └─ No ↓
  │
  └─ Is it INLINE TEXT (pasted content, agent output, conversation excerpt)?
      └─ Yes → Apply fidelity rules directly
               Use extractive method: identify key passages first
               Produce structured output
```

## Delegation Decision

When the summarization task is autonomous (user asks to summarize something and move on), delegate to a specialized agent:

| Source Type | Agent | When to Use |
|-------------|-------|-------------|
| File(s) | @file-summarizer | Summarizing files without immediate follow-up questions |
| URL(s) | @url-summarizer | Summarizing web content autonomously |
| Image(s) | @image-summarizer | Describing visual content autonomously |

When the summarization is part of the current conversation flow (user wants to discuss the content), apply the relevant skill methodology directly rather than delegating.

**Orchestrator relay**: When receiving results from any summarizer agent, the orchestrator MUST follow the agent-result-relay skill to preserve counts, failure reasons, and structured output. Do not re-summarize agent summaries.

## Team Coordination (Multi-Source)

When the Teammate tool is available and 3+ sources require summarization, the model SHOULD use agent teams instead of sequential subagents. If the Teammate tool is not available, fall back to subagent delegation.

### Workflow

1. **Create team** - `Teammate({ operation: "spawnTeam", team_name: "summarize-{task-id}" })`
2. **Create tasks** - One TaskCreate per source, all independent (no dependencies)
3. **Spawn teammates** - One per source, using the appropriate agent type:

```text
Task({
  team_name: "summarize-{task-id}",
  name: "source-1",
  subagent_type: "file-summarizer",  // or url-summarizer, image-summarizer
  prompt: "Summarize [source path]. Format: {format_id}. When done, send findings to team-lead via Teammate write. If you find information that contradicts another source, message that teammate directly.",
  run_in_background: true
})
```

4. **Collect results** - Leader receives findings via inbox messages
5. **Synthesize** - Leader applies multi-source-synthesis skill to the collected findings
6. **Cleanup** - Request shutdown for all teammates, then cleanup

### When NOT to Use Teams

- Fewer than 3 sources (subagent overhead is lower)
- Sources have sequential dependencies (one source references another)
- User is in a conversational flow and wants to discuss each source

### Fidelity Rules Still Apply

All fidelity rules apply identically to teammate output. The SubagentStop hook validates teammate summaries the same way it validates subagent summaries.

## Fidelity Rules (Mandatory)

These rules apply to ALL summarization regardless of source type. See [Fidelity Rules](./references/fidelity-rules.md) for full details.

**Rule 1: Read Before Summarizing** - Read actual content. Never guess from filenames or paths.

**Rule 2: Extract Before Abstracting** - Pull quotes/passages first, then summarize from extracts.

**Rule 3: Preserve Counts and Specifics** - Keep exact numbers. "7 of 10" not "most."

**Rule 4: Distinguish Absence from Nonexistence** - "Not mentioned in source" not "doesn't exist."

**Rule 5: No Lossy Re-Summarization** - When relaying agent results, relay counts and references. Do not summarize the summary.

**Rule 6: State Confidence Explicitly** - Every summary includes confidence level with rationale.

**Rule 7: Structured Output Always** - Use the format defined in [Structured Summary](./templates/structured.md).

## Output Format

All summaries MUST use structured markdown with YAML frontmatter. See [Structured Summary](./templates/structured.md) for the complete specification.

Required sections in every summary:

1. **YAML frontmatter** - source_type, source_path, method, confidence, word counts
2. **Summary** - the condensed content (BLUF style)
3. **What Was Found** - items discovered with source references
4. **What Was NOT Found** - items searched for but absent
5. **Uncertain** - ambiguous items requiring interpretation
6. **Sources** - full attribution with access dates

## Anti-Patterns

The model MUST NOT:

- Summarize a file based on its name without reading it
- Summarize a URL based on its domain without fetching it
- Describe an image based on its filename without viewing it
- Re-summarize a sub-agent's summary (relay instead)
- Upgrade "not found" to "doesn't exist"
- Drop counts ("7 of 10" → "most")
- Omit the "What Was NOT Found" section
- Present uncertain information as definitive
- Summarize from excerpts (head/tail/grep) without disclosing the limitation
