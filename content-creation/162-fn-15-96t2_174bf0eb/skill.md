# fn-15-96t.2 Create plan-sync agent

## Description
TBD

## Acceptance
- [ ] TBD

## Done summary
Created plan-sync agent (plugins/flow-next/agents/plan-sync.md) with phased workflow: re-anchor on completed task, explore implementation, identify drift, check downstream tasks, update affected specs. Uses opus model with Read/Grep/Glob/Edit tools, disallows Bash/Write/Task.
## Evidence
- Commits: f3ee708efc5bfd2056c6780b78511e8bd63c45bf
- Tests:
- PRs: