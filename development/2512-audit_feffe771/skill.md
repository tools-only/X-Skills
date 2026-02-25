# Audit Workflow

Loaded by: `/rwr:audit` command
Orchestrator: Claude (reads this workflow and executes steps)

## Step 1 — Classify Task

Determine task type from $ARGUMENTS:

```mermaid
flowchart TD
    Start([$ARGUMENTS]) --> Q{Signal words?}
    Q -->|"drift, out of date, stale docs, undocumented, documented but, doesn't match code, audit"| Drift[Task type: drift-audit]
    Q -->|"sync, update docs, docs behind, after refactor, after implementing"| Sync[Task type: documentation-sync]
    Q -->|"freshness, stale headers, governance, track changes, last updated"| Fresh[Task type: freshness]
    Q -->|Ambiguous| Ask[Ask user — see disambiguation below]
```

Disambiguation prompt when ambiguous:

> "Is this (a) checking if docs match code, (b) updating docs after code changed, or (c) tracking doc freshness?"

## Step 2 — Read Agent Protocol

Before spawning, read the target agent's file to understand its exact required inputs and output format:

- drift-audit: Read `plugins/development-harness/agents/doc-drift-auditor.md`
- documentation-sync: Read `plugins/development-harness/agents/service-docs-maintainer.md`
- freshness: Read `/home/ubuntulinuxqa2/.claude/agents/doc-freshness-guardian.md`

## Step 3 — Spawn Agent

Construct prompt using exact input format the agent expects (from Step 2 read).

For drift-audit:

```text
Task(
  subagent_type="development-harness:doc-drift-auditor",
  prompt="<task description>

Scope: <file paths or directories>
Project root: <path>"
)
```

For documentation-sync:

```text
Task(
  subagent_type="development-harness:service-docs-maintainer",
  prompt="<description of what code changed>

Affected files: <list>"
)
```

Note: service-docs-maintainer does NOT write a summary file — its output is response text only.

For freshness:

```text
Task(
  subagent_type="doc-freshness-guardian",
  prompt="<task description>

Files to check: <list>"
)
```

## Step 4 — Handle Return

```mermaid
flowchart TD
    Return([Agent returns]) --> Q{Agent type?}
    Q -->|drift-audit STATUS: DONE| Report[Report findings to user]
    Report --> Chain{Markdown files modified?}
    Chain -->|Yes| Validate[Chain to validate workflow\nLoad: plugins/the-rewrite-room/the-rewrite-room/workflows/validate.md]
    Chain -->|No| Done([Done])
    Q -->|drift-audit STATUS: BLOCKED| ShowNeeded[Show NEEDED list to user]
    ShowNeeded --> Ask[Ask user to provide missing inputs]
    Ask --> Retry[Return to Step 3 with complete inputs]
    Q -->|documentation-sync completes| ShowSummary[Show summary of docs changed]
    ShowSummary --> Offer["Offer: Run drift audit to verify docs now match code? (y/n)"]
    Q -->|freshness returns| ShowFreshness[Show freshness report]
    ShowFreshness --> Highlight[Highlight red greater than 90d stale items first]
```

## Output Contract

```text
STATUS: DONE|BLOCKED|FAILED
SUMMARY: [1-2 sentences — what was audited/synced, key findings count]
ARTIFACTS: [files created/modified with relative paths, or "none"]
VALIDATION:
  - citation-check: PASS|FAIL (drift-audit only — all findings have file:line evidence)
  - link-check: PASS|FAIL (if markdown files modified)
NOTES: [only if needed]
```
