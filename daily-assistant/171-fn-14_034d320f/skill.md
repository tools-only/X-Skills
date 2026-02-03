# fn-1.4 Memory-scout subagent

## Description
Create memory-scout subagent for retrieving relevant memory during planning and work phases. One agent, two call sites:

1. `/flow-next:plan` Step 1 - parallel with repo-scout, practice-scout, docs-scout
2. `/flow-next:work` Phase 1 - during re-anchor after flowctl show/cat

Agent reads `.flow/memory/` files and returns only relevant entries (not everything).

Call pattern (matches existing scouts):
```
Task flow-next:memory-scout(<request>)
```

## Acceptance
- [ ] `agents/memory-scout.md` created with agent definition
- [ ] Agent reads pitfalls.md, conventions.md, decisions.md
- [ ] Returns only semantically relevant entries
- [ ] Added to plan skill Step 1 (gated by memory.enabled)
- [ ] Added to work skill Phase 1 re-anchor (gated by memory.enabled)
- [ ] Works with empty memory (returns nothing gracefully)
- [ ] Handles large memory files (selective return)

## Done summary
- Created agents/memory-scout.md with haiku model
- Updated plan skill steps.md to call memory-scout in Step 1
- Updated work skill phases.md to call memory-scout in Phase 3
- All calls gated by `$FLOWCTL config get memory.enabled`

Why:
- Single agent for both plan and work phases
- Uses haiku for efficiency (simple lookup task)
- Returns only relevant entries, not entire memory

Verification:
- Agent file follows existing scout format
- flowctl validate --all passes
## Evidence
- Commits: 73fa01c602619ea30031cc17e7098c07c3db23c5
- Tests: flowctl validate --all
- PRs: