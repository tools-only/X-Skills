---
description: "Groom backlog items — trigger /groom-backlog-item <title|section|all> — runs RT-ICA per item then spawns @backlog-item-groomer agents to discover research, skills, agents, prior work, and dependencies. Produces context manifests and grooming report. Use when preparing backlog items for planning or execution."
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

### Step 3: RT-ICA Assessment Per Item

Perform Reverse Thinking — Information Completeness Assessment before spawning groomer agents. This directs the groomer's discovery toward filling gaps rather than broad search.

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

- **AVAILABLE**: Explicitly stated in item description or research questions
- **DERIVABLE**: Safely inferable from codebase context (state basis)
- **MISSING**: Not present, not safely inferable — groomer must find it

Pass the RT-ICA summary to the groomer alongside item details.

### Step 4: Spawn Groomer Agents

**Single item** — invoke `@backlog-item-groomer` directly, passing item details and RT-ICA summary.

**Multiple items** — spawn parallel Task agents (max 5 concurrent; batch in waves if more):

```text
Task(
  subagent_type: "general-purpose",
  prompt: "Act as @backlog-item-groomer. Groom this item:\n{item details}\n\nRT-ICA Assessment:\n{rt-ica summary}",
  model: "haiku"
)
```

### Step 5: Collect and Report

Gather context manifests. Produce grooming report:

```markdown
# Backlog Grooming Report

**Date**: {YYYY-MM-DD}
**Items groomed**: {count}
**Arguments**: {original arguments}

## Summary

| Item | RT-ICA | Research Found | Skills | Agents | Blockers |
|------|--------|----------------|--------|--------|----------|
| {title} | {APPROVED/BLOCKED} | {count} | {count} | {count} | {count} |

## Individual Manifests

### {Item title}
{manifest from agent}

## RT-ICA Results

### BLOCKED Items
- {item title}: {list of missing inputs}

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

- RT-ICA summary included for each item
- Groomer agent(s) received RT-ICA context
- Report contains RT-ICA Results section
- Cross-item findings present (if multiple items groomed)
