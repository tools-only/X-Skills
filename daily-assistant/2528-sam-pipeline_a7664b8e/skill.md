# SAM 7-Stage Pipeline

Canonical pipeline for feature development. Each stage produces a file-based artifact; no stage relies on conversation memory.

---

## Pipeline Overview

```mermaid
flowchart TD
    Start([Feature Request]) --> S1[Stage 1 — Discovery]
    S1 -->|ARTIFACT:DISCOVERY| HT1{Human Touchpoint?}
    HT1 -->|ARL — unbound constraints| Human1[Escalate]
    HT1 -->|ARL — bound constraints| S2[Stage 2 — Planning + RT-ICA]
    Human1 --> S2
    S2 -->|ARTIFACT:PLAN| S3[Stage 3 — Context Integration]
    S3 --> S4[Stage 4 — Task Decomposition]
    S4 -->|ARTIFACT:TASK| HT2{Human Touchpoint?}
    HT2 -->|ARL — high complexity| Human2[Escalate]
    HT2 -->|ARL — routine| S5[Stage 5 — Execution]
    Human2 --> S5
    S5 -->|ARTIFACT:EXECUTION| S6[Stage 6 — Forensic Review]
    S6 -->|NEEDS_WORK| S5
    S6 -->|COMPLETE| S7[Stage 7 — Final Verification]
    S7 -->|NOT_CERTIFIED| S4
    S7 -->|CERTIFIED| Done([Feature Complete])
```

---

## Stages

| Stage | Name | Input | Output |
|-------|------|-------|--------|
| S1 | Discovery | User request, problem statement | ARTIFACT:DISCOVERY(SCOPE:...) |
| S2 | Planning + RT-ICA | Discovery artifacts | ARTIFACT:PLAN(SCOPE:...) |
| S3 | Context Integration | Plan + codebase | ARTIFACT:PLAN contextualized |
| S4 | Task Decomposition | Contextualized plan | ARTIFACT:TASK(TASK:...) per task |
| S5 | Execution | Single task file | ARTIFACT:EXECUTION(TASK:...) |
| S6 | Forensic Review | Execution + task + plan | ARTIFACT:REVIEW — COMPLETE / NEEDS_WORK |
| S7 | Final Verification | All completed tasks + goals | ARTIFACT:VERIFICATION — CERTIFIED / NOT_CERTIFIED |

---

## Loop Limits

- **NEEDS_WORK**: 3 iterations per task before human escalation
- **NOT_CERTIFIED**: 2 iterations before human escalation

---

## Artifact Flow

```text
User Request → DISCOVERY → PLAN → PLAN (contextualized) → TASK(s) → EXECUTION(s) → REVIEW(s) → VERIFICATION
```

---

## Source

- [sam-definition.md](../../.claude/skills/work-backlog-item/references/sam-definition.md)
- [default-development-flow.md](../../plugins/development-harness/skills/development-harness/references/default-development-flow.md)
