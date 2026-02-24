---
name: groom-backlog-item
description: Groom backlog items — trigger /groom-backlog-item <title|section|all> — fact-checks item claims against primary sources, runs RT-ICA per item, then spawns @backlog-item-groomer agents. Writes groomed content into per-item files in .claude/backlog/. Use when preparing backlog items for planning or execution.
argument-hint: <item-title-or-section-or-all>
user-invocable: true
---
# Groom Backlog Item

Orchestrate backlog grooming: parse arguments, assess information completeness via RT-ICA, spawn discovery agents, write groomed content into per-item files.

## Arguments

`$ARGUMENTS` accepts:

- **Title substring** — e.g., `Error Recovery` — grooms matching item (case-insensitive)
- **Section** — `P0`, `P1`, `P2`, or `Ideas` — grooms all items in that section
- **`all`** — grooms all items across P0, P1, P2, Ideas (parallel agents)

## Workflow

### Step 1: Parse Arguments and Load Backlog

Read `.claude/BACKLOG.md`. Identify target items based on argument type above.

### Step 2: Validity Check (Pre-Groom Gate)

Before fact-checking or grooming, verify each item is still valid work:

1. **Is the job still valid?** — Scope, priority, or context may have changed. Ask or infer: does this item still belong in the backlog?
2. **Can the problem be replicated?** — If the item describes a bug or fix, confirm the issue still exists. If it cannot be reproduced, consider resolving or closing.
3. **Is this local file stale?** — If the item has a GitHub issue (`metadata.issue` or index link `#N`), fetch the issue state via `gh issue view N --json state` (use `-R Jamie-BitFlight/claude_skills` if in a proxy environment). If the issue is **closed**, the local file is a stale remnant of work already done. Do **not** groom. Instead:
   - Recommend: `backlog close "{title}" --plan <path> --checklist-pass --cleanup` (if completed) or `backlog resolve "{title}" --reason "..." --cleanup` (if obsolete)
   - Skip grooming for that item; move to the next

If any check fails, skip grooming for that item and report. Only proceed to Step 3 for items that pass.

### Step 3: Extract Item Details

For each target item, extract: title, description, research-first questions (if present), source, suggested location.

### Step 4: Fact-Check Item Claims

Invoke the `fact-check` skill on each target item to verify factual claims against primary sources **before** running RT-ICA or spawning groomer agents. This prevents unverified or refuted assertions from entering the planning context.

```text
Skill(command: "fact-check", args: "{item title}")
```

The `fact-check` skill spawns `@fact-checker` agents that MUST retrieve evidence via `WebFetch`, `WebSearch`, or `gh`. Training data recall is not accepted as evidence.

After each run, collect the verdict summary:

```text
Fact-Check Summary: {item title}
Claims checked: {N}
VERIFIED: {N} | REFUTED: {N} | INCONCLUSIVE: {N}
Refuted claims:      [{list of claim texts — each becomes a MISSING condition in Step 4}]
Inconclusive claims: [{list of claim texts — flag as unverified DERIVABLE in Step 4}]
Citations:           [{VERIFIED claims cite their primary sources}]
```

**Multiple items** — invoke `fact-check` for each item sequentially (respect the wave-of-5 concurrency limit inside `fact-check` itself). Do not batch items into a single `fact-check` call.

Pass the fact-check summary forward to Step 5.

### Step 5: RT-ICA Assessment Per Item

Perform Reverse Thinking — Information Completeness Assessment using both the item details **and** the fact-check verdicts from Step 4. This directs the groomer's discovery toward filling gaps rather than broad search.

For each item, produce:

```text
RT-ICA: {item title}
Goal: {one sentence — what completing this item achieves}
Conditions:
1. {condition} | Status: {AVAILABLE|DERIVABLE|MISSING} | Info needed: {what}
...
Decision: {APPROVED|BLOCKED}
Missing: {list of missing inputs, or "None"}
```

- **AVAILABLE**: Explicitly stated in item description or research questions AND fact-check verdict (Step 4) is VERIFIED or not applicable
- **DERIVABLE**: Safely inferable from codebase context (state basis); fact-check verdict is INCONCLUSIVE
- **MISSING**: Not present, not safely inferable — OR fact-check verdict is REFUTED (the stated condition is false and the correct state is unknown)

REFUTED claims from Step 4 MUST be listed as MISSING conditions. A REFUTED claim is not a valid basis for any AVAILABLE or DERIVABLE status.

Pass the RT-ICA summary and fact-check summary to the groomer alongside item details.

**ARL human-probing integration:** When RT-ICA returns BLOCKED or MISSING conditions, the context manifest can include `invisible_knowledge_prompts` — questions to ask the human before planning (e.g., "What went wrong in the past?", "What references are essential?"). See [.claude/docs/sdlc-layers/arl-human-probing-design.md](../../docs/sdlc-layers/arl-human-probing-design.md).

### Step 6: Spawn Groomer Agents

**Single item** — invoke `@backlog-item-groomer` directly, passing item details, RT-ICA summary, and fact-check summary.

**Multiple items** — spawn parallel Task agents (max 5 concurrent; batch in waves if more):

```text
Task(
  subagent_type: "general-purpose",
  prompt: "Act as @backlog-item-groomer. Groom this item and output groomed content in the standard template format (see .claude/docs/backlog-item-groomed-schema.md). Output only the groomed body (no ## Groomed header).\n\nItem:\n{item details}\n\nRT-ICA Assessment:\n{rt-ica summary}\n\nFact-Check Verdicts:\n{fact-check summary}",
  model: "haiku"
)
```

### Step 7: Write Groomed Content to Item Files

For each item, write groomed content into the per-item file via the backlog script. Prefer incremental updates so sections (Fact-Check, RT-ICA, groomed subsections) are appended as they become available. GitHub is canonical: when the item has an issue, the backlog script syncs groomed content to the GitHub issue body.

**Preferred: incremental section updates**

After each step, call the backlog script with `--section` and `--content`:

```text
# After Step 3 (fact-check)
backlog groom "{item title}" --section "Fact-Check" --content "{fact-check summary}"

# After Step 4 (RT-ICA)
backlog groom "{item title}" --section "RT-ICA" --content "{rt-ica summary}"

# After Step 5 (groomer output) — full groomed body or subsections
backlog groom "{item title}" --section "Reproducibility" --content "{reproducibility section}"
# ... or for full groomed body:
backlog groom "{item title}" --groomed-content "{full groomed body}"
```

**Alternative: full content**

```text
backlog groom "{item title}" --groomed-content "{full groomed body}"
# Or from file:
backlog groom "{item title}" --groomed-file {path}
# Or from stdin:
backlog groom "{item title}" < {groomed_file}
```

**Valid section names** — top-level: `Fact-Check`, `RT-ICA`. Groomed subsections: `Reproducibility`, `Priority`, `Impact`, `Scope`, `Output / Evidence`, `Dependencies`, `Research`, `Skills`, `Agents`, `Prior Work`, `Files`, `Decision`.

The backlog script updates `.claude/backlog/{priority}-{slug}.md` with merged sections, sets `groomed` in frontmatter, and syncs to the GitHub issue when the item has one.

**Bulk grooming (multiple items)** — when grooming 2+ items, optionally persist a session summary to `.claude/grooming-sessions/{YYYY-MM-DD}.md`:

```markdown
# Grooming Session {YYYY-MM-DD}

**Items groomed**: {count}
**Arguments**: {original arguments}

## Summary

| Item | Fact-Check | RT-ICA | Written |
|------|------------|--------|---------|
| {title} | {V}/{R}/{I} | {APPROVED/BLOCKED} | ✓ |

## Cross-Item Findings

### Shared Dependencies
- {items multiple backlog items depend on}

### Suggested Groupings
- {items that could be worked together}

### Research Gaps
- {topics needing research}
```

Per-item groomed content lives in each item file; this session file holds only metadata and cross-item findings.

## Example Invocations

```text
/groom-backlog-item Error Recovery
/groom-backlog-item P1
/groom-backlog-item all
```

## Completion Criteria

- Validity check (job still valid, problem reproducible, local file not stale) before grooming
- Fact-check run for each item before RT-ICA (training data not used as evidence)
- Fact-check verdicts passed into RT-ICA conditions (REFUTED → MISSING)
- RT-ICA summary included for each item
- Groomer agent(s) received RT-ICA context and fact-check verdicts
- Groomed content written via `backlog groom` (prefer `--section`/`--content` incremental updates; `--groomed-content` or stdin for full body)
- When item has GitHub issue, groomed content synced to issue body
- Bulk session summary optionally saved to `.claude/grooming-sessions/{date}.md` when grooming multiple items
