---
title: "Token Pruning Efficiency Improvements – Plan"
phase: Plan
date: "2026-02-02T00:00:00Z"
owner: "claude"
parent_research: "memory-bank/research/2026-02-02_token_estimation_overhead.md"
git_commit_at_plan: ""
tags: [plan, performance, pruning]
---

## Goal

Reduce pruning overhead in `src/tunacode/core/agents/resume/prune.py` by eliminating redundant token estimation, while preserving current pruning behavior.

**Non-goals:**
- Changing pruning thresholds or user-visible behavior
- Introducing new dependencies or schema changes without approval
- Broad refactors outside `prune.py` and token estimation helpers

## Scope & Assumptions

**In scope:**
- Reuse precomputed token counts for pruning mutations
- Precompute placeholder token counts once
- Stop token estimation in Phase 1 once thresholds are satisfied (continue scanning to find parts)
- Optional: lightweight token caching strategy (if safe for message schema)

**Out of scope:**
- Replacing `estimate_tokens()` heuristic
- Changes to message history structure unless explicitly approved
- UI or telemetry changes

**Assumptions:**
- Token estimation remains model-agnostic
- Existing pruning thresholds remain valid
- Current tests are minimal; manual validation will be needed

## Deliverables (DoD)

1. **No redundant token estimation** in Phase 4 (reuses Phase 1 values)
2. **Placeholder token count cached** as a constant
3. **Stop token estimation** in Phase 1 after thresholds are satisfied while still scanning for parts to prune
4. **No new mypy errors** (baseline maintained)
5. **Ruff clean** (`uv run ruff check .`)

## Readiness (DoR)

- [x] Research doc complete
- [ ] Target function behavior validated with sample history
- [ ] Agreement on token caching approach (if any)

## Milestones

- **M1:** Remove redundant token estimation (Phase 4 reuse)
- **M2:** Precompute placeholder token count constant
- **M3:** Stop token estimation after threshold while continuing Phase 1 scan
- **M4:** Optional caching strategy (decision gated)
- **M5:** Validation (manual + targeted tests if available)

## Work Breakdown (Tasks)

| Task | ID | Summary | File | Acceptance Tests |
|------|----|---------|------|------------------|
| T1 | tk-new | Reuse stored token counts during pruning mutation | `src/tunacode/core/agents/resume/prune.py` | Manual diff of pruned output (before/after) |
| T2 | tk-new | Add `PRUNE_PLACEHOLDER_TOKENS` constant | `src/tunacode/core/agents/resume/prune.py` | Placeholder count unchanged |
| T3 | tk-new | Stop Phase 1 token estimation after `protect + minimum` while still scanning for parts | `src/tunacode/core/agents/resume/prune.py` | Same pruning boundary on sample data |
| T4 | tk-new | Decide on caching: add per-part token cache or skip | `prune.py` or `token_counter.py` | Explicit decision captured |

## Risks & Mitigations

| Risk | Impact | Likelihood | Mitigation | Trigger |
|------|--------|------------|------------|---------|
| Skipping scan would change prune boundary | Medium | Low | Ensure scan continues; validate with fixed histories | Different prune count |
| Token cache mutates message schema | Medium | Medium | Prefer side map keyed by id; gate with approval | Reviewer concern |
| Placeholder tokens drift | Low | Low | Compute once with existing estimator | Test mismatch |

## Test Strategy

- **Manual validation:**
  - Run pruning on a sample history before/after and compare:
    - total reclaimed tokens
    - pruned parts count
  - Ensure pruned text equals `PRUNE_PLACEHOLDER`
- **Automated (if time):**
  - Add targeted unit test for `prune_old_tool_outputs()` with fixed content
- **Lint:** `uv run ruff check .`

## References

- Research doc: `memory-bank/research/2026-02-02_token_estimation_overhead.md`
- Target file: `src/tunacode/core/agents/resume/prune.py`
- Token estimator: `src/tunacode/utils/messaging/token_counter.py`

## Tickets Created (max 5)

| Ticket ID | Title | Priority | Status |
|-----------|-------|----------|--------|
| tun-token-1 | Reuse precomputed token counts in pruning | 1 | open |
| tun-token-2 | Add placeholder token constant | 2 | open |
| tun-token-3 | Early termination in tool part collection | 2 | open |
| tun-token-4 | Decide on token caching strategy | 3 | open |
