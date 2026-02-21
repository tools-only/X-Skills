---
adapter: contextual-ai-documentation-optimizer-adapter
wraps: plugins/plugin-creator/agents/contextual-ai-documentation-optimizer.md
normalizes_to: optimization-block-v1
version: "1.0"
---

# contextual-ai-documentation-optimizer Adapter

## Wrapped Component

`plugins/plugin-creator/agents/contextual-ai-documentation-optimizer.md`

This agent produces a multi-section output: RT-ICA assessment, optimized content, CoVe verification, and token impact report. Output is a STATUS: DONE or BLOCKED block with the optimized content inline or written to file.

## Native Output Contract

The agent produces:

1. RT-ICA pre-check assessment
2. Optimized file content (written in-place or to specified path)
3. CoVe post-check verification (did the output meet the objectives?)
4. Token impact report (before/after token counts)
5. STATUS: DONE or BLOCKED

## Adapter Transformation

The adapter maps the native output to `optimization-block-v1`:

```text
STATUS: DONE|BLOCKED|FAILED
SUMMARY: [what was optimized — extracted from agent's STATUS summary]
ARTIFACTS:
  - {optimized-file-path}
VALIDATION:
  - frontmatter-validator: PASS|FAIL
  - prompt-structure-validator: PASS|FAIL
DIFF:
  tokens_before: N
  tokens_after: N
  delta: +N / -N
  prohibitions_converted: N
NOTES: [CoVe verification notes, if any issues found]
```

## Validators Added by Adapter

- `frontmatter-validator` — re-runs after optimization to confirm frontmatter is still valid
- `prompt-structure-validator` — checks that prohibition patterns were converted to positive alternatives

## BLOCKED Mapping

If the agent returns BLOCKED (typically due to RT-ICA pre-check failure), the adapter preserves the BLOCKED status. The RT-ICA failure reason is placed in NOTES. The user must provide the missing input before re-running.

## Token Impact Data

The adapter extracts token counts from the agent's token impact report. If the optimized content increased token count by more than 10%, the NOTES field includes a flag: "Token count increased — review before committing."
