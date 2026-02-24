# Artifact Conventions

SAM artifact naming, file layout, and cross-referencing. Artifacts are the only communication channel between stages.

---

## Token Pattern

**Format:** `ARTIFACT:{TYPE}({SCOPE_OR_ID})`

**Types:**

- DISCOVERY, PLAN, CONTEXT, TASK, EXECUTION, REVIEW, VERIFICATION

---

## File Layout

All artifacts in `.planning/harness/`:

```text
.planning/harness/
  discovery-{feature-slug}.md
  plan-{feature-slug}.md
  context-{feature-slug}.md
  task-{task-id}-{task-slug}.md
  execution-{task-id}-{task-slug}.md
  review-{feature-slug}.md
  verification-{feature-slug}.md
```

---

## Coexistence

- `.planning/gsd/` — GSD planning
- `.planning/backlog/` — Backlog items
- `.planning/` root — Shared planning docs

---

## Source

- [artifact-conventions.md](../../plugins/development-harness/skills/development-harness/references/artifact-conventions.md)
