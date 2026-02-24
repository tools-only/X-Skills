---
name: groom-backlog-item
description: Groom backlog items — trigger /groom-backlog-item <title|section|all> — fact-checks item claims against primary sources, runs RT-ICA per item, then spawns @backlog-item-groomer agents to discover research, skills, agents, prior work, and dependencies. Produces context manifests and grooming report. Use when preparing backlog items for planning or execution.
argument-hint: <item-title-or-section-or-all>
user-invocable: true
---
# Groom Backlog Item

Orchestrate backlog grooming: parse arguments, assess information completeness via RT-ICA, spawn discovery agents, produce report.

## Arguments

`$ARGUMENTS` accepts:

- **Title substring** — e.g., `Error Recovery` — grooms matching item (case-insensitive)
- **Section** — `P0`, `P1`, `P2`, or `Ideas` — grooms all items in that section
- **`all`** — grooms all items across P0, P1, P2, Ideas (parallel agents)

## Workflow

### Step 1: Parse Arguments and Load Backlog

Read `.claude/BACKLOG.md`. Identify target items based on argument type above.

### Step 2: Extract Item Details

For each target item, extract: title, description, research-first questions (if present), source, suggested location.

### Step 3: Fact-Check Item Claims

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

Pass the fact-check summary forward to Step 4.

### Step 4: RT-ICA Assessment Per Item

Perform Reverse Thinking — Information Completeness Assessment using both the item details **and** the fact-check verdicts from Step 3. This directs the groomer's discovery toward filling gaps rather than broad search.

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

- **AVAILABLE**: Explicitly stated in item description or research questions AND fact-check verdict is VERIFIED or not applicable
- **DERIVABLE**: Safely inferable from codebase context (state basis); fact-check verdict is INCONCLUSIVE
- **MISSING**: Not present, not safely inferable — OR fact-check verdict is REFUTED (the stated condition is false and the correct state is unknown)

REFUTED claims from Step 3 MUST be listed as MISSING conditions. A REFUTED claim is not a valid basis for any AVAILABLE or DERIVABLE status.

Pass the RT-ICA summary and fact-check summary to the groomer alongside item details.

**ARL human-probing integration:** When RT-ICA returns BLOCKED or MISSING conditions, the context manifest can include `invisible_knowledge_prompts` — questions to ask the human before planning (e.g., "What went wrong in the past?", "What references are essential?"). See [.claude/docs/sdlc-layers/arl-human-probing-design.md](../docs/sdlc-layers/arl-human-probing-design.md).

### Step 5: Spawn Groomer Agents

**Single item** — invoke `@backlog-item-groomer` directly, passing item details, RT-ICA summary, and fact-check summary.

**Multiple items** — spawn parallel Task agents (max 5 concurrent; batch in waves if more):

```text
Task(
  subagent_type: "general-purpose",
  prompt: "Act as @backlog-item-groomer. Groom this item:\n{item details}\n\nRT-ICA Assessment:\n{rt-ica summary}\n\nFact-Check Verdicts:\n{fact-check summary}",
  model: "haiku"
)
```

### Step 6: Collect and Report

Gather context manifests. Produce grooming report:

```markdown
# Backlog Grooming Report

**Date**: {YYYY-MM-DD}
**Items groomed**: {count}
**Arguments**: {original arguments}

## Summary

| Item | Fact-Check | RT-ICA | Research Found | Skills | Agents | Blockers |
|------|------------|--------|----------------|--------|--------|----------|
| {title} | {V}/{R}/{I} | {APPROVED/BLOCKED} | {count} | {count} | {count} | {count} |

## Individual Manifests

### {Item title}
{manifest from agent}

## Fact-Check Results

### Refuted Claims
- {item title}: {claim text} — REFUTED by {source URL}

### Inconclusive Claims
- {item title}: {claim text} — INCONCLUSIVE: {what additional verification is needed}

### Verified Claims
- {item title}: {count} claims verified against primary sources

## RT-ICA Results

### BLOCKED Items
- {item title}: {list of missing inputs, including any from refuted claims}

### APPROVED Items
- {item title}: {count} conditions verified

## Cross-Item Findings

### Shared Dependencies
- {items multiple backlog items depend on}

### Suggested Groupings
- {items that could be worked together}

### Research Gaps
- {topics needing research — methodology: [stateless-agent-methodology](https://github.com/bitflight-devops/stateless-agent-methodology)}
```

If grooming multiple items, offer to save report to `.claude/grooming-reports/grooming-{YYYY-MM-DD}.md`.

## Example Invocations

```text
/groom-backlog-item Error Recovery
/groom-backlog-item P1
/groom-backlog-item all
```

## Completion Criteria

- Fact-check run for each item before RT-ICA (training data not used as evidence)
- Fact-check verdicts passed into RT-ICA conditions (REFUTED → MISSING)
- RT-ICA summary included for each item
- Groomer agent(s) received RT-ICA context and fact-check verdicts
- Report contains Fact-Check Results section and RT-ICA Results section
- Cross-item findings present (if multiple items groomed)
