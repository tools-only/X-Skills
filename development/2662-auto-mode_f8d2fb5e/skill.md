# Auto Mode Rules

When `$0` is `--auto`, the following substitutions apply at every interactive decision point:

| Normal behaviour | `--auto` substitution |
|---|---|
| No title given (`$1` is empty) | Scan `.claude/backlog/` per-item files for P0 then P1 items with status `needs-grooming` or `groomed`. Log `[AUTO] No title — auto-selected: {title}` and proceed. If none found, log `[AUTO] STOP — no open P0/P1 items found` and stop. |
| Step 1b: issue not found | Log `[AUTO] STOP — Issue #N not found`, stop |
| Step 1: zero matches → ask user to create | Auto-invoke `create-backlog-item --auto {title}`, log `[AUTO] No item found — invoking create-backlog-item --auto` |
| Step 1: multiple matches → ask user to pick | Log `[AUTO] Multiple matches — picking first: {title}`, proceed with first match |
| Step 2.5: offer GitHub issue for P0/P1 | Log `[AUTO] Skipping GitHub issue offer`, continue without issue |
| Step 2.5: ask milestone assignment | Log `[AUTO] Skipping milestone assignment`, skip |
| RT-ICA BLOCKED | Log `[AUTO] STOP — RT-ICA BLOCKED: {missing inputs}`, stop (cannot resolve without human) |
| Any other `AskUserQuestion` | Log `[AUTO] Decision: {chosen option} — reason: {evidence}`, proceed with logged choice |

`--auto` does NOT change the behaviour of Steps 3–8 (grooming, RT-ICA evaluation, SAM planning, backlog update) — those are already agent-executable without human input.
