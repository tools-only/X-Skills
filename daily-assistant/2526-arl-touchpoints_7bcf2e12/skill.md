# Human Touchpoint Model

ARL-derived escalation decision model. Not every stage transition requires human approval. The harness uses constraint analysis to escalate only when the agent cannot safely proceed.

---

## Escalation Format

When escalating, present:

1. **What stage** produced the escalation
2. **What triggered** it (specific constraint, risk factor, or loop limit)
3. **What the agent knows** (bound constraints, available context)
4. **What the agent does not know** (unbound constraints, missing information)
5. **Decision options** — concrete choices the human can make to unblock

**Never** present vague "please review" requests.

---

## Pre-Scheduled Gates

**Gate 1 — After S1 Discovery, before S2 Planning:**

- Triggered when: unbound constraints, domain knowledge gaps
- Skipped when: all constraints bound, sufficient context

**Gate 2 — After S4 Task Decomposition, before S5 Execution:**

- Triggered when: high complexity, novel architecture
- Skipped when: routine patterns, existing codebase precedent

---

## Dynamic Escalation Points

- **NEEDS_WORK loop limit** — 3 iterations in S6 without resolution
- **NOT_CERTIFIED loop limit** — 2 iterations in S7 without resolution
- Agent failure, quality gate cascade failure, contradiction detected

---

## Constraint Types

| Type | Meaning | Action |
|------|---------|--------|
| **Bound** | All information available, deterministic | Proceed autonomously |
| **Unbound** | Information missing, ambiguous, or requires external knowledge | Escalate |
| **Mixed** | Some bound, some unbound | Isolate unbound; proceed with bound if possible |

---

## Source

- [human-touchpoint-model.md](../../plugins/development-harness/skills/development-harness/references/human-touchpoint-model.md)
