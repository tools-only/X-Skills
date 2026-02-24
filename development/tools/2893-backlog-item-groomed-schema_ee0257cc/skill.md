# Backlog Item Groomed Schema

**Purpose**: Define the structure of groomed content written into backlog item files. Grooming transforms items from vague ("this problem happens") to ready ("reproducible, impact known, desired outcome clear").

**Location**: Groomed content lives in the **body** of `.claude/backlog/{priority}-{slug}.md`. Frontmatter uses the research-style `metadata:` block (aligned with `./research/` entries). Body has no duplication of frontmatter — only extra fields when present, plus `## Groomed` when groomed.

---

## Frontmatter (Research-Style)

Aligned with [research/README.md](../../research/README.md) metadata format:

```yaml
---
name: "Item title"
description: "Main problem statement or description"
metadata:
  topic: "slug-from-title"
  source: "Source description"
  added: "YYYY-MM-DD"
  priority: P0|P1|P2|Ideas
  type: Feature|Bug|Refactor|Docs|Chore
  status: open|in-progress|done|resolved
  # optional:
  issue: "#N"
  plan: "path/to/plan.md"
  groomed: "YYYY-MM-DD"
  layer: "0"|"1"|"2"
  language: python|typescript|...
  stack: fastapi|python-cli|...
---
```

---

## Body (No Duplication)

Body is **empty** when un-groomed and no extra details. When present, body contains only:

- **Suggested location** — when applicable
- **Research first** — when applicable
- **Decision needed** — when applicable
- **Files** — when applicable
- **## Groomed (YYYY-MM-DD)** — when groomed (see below)

Do **not** duplicate `name`, `description`, `source`, `added`, `priority`, `type`, `issue` in the body — they live in frontmatter.

---

## Groomed Sections (Body, under ## Groomed)

| Section | Purpose | Required |
|---------|---------|----------|
| **Reproducibility** | Steps x, y, z to replicate the problem | When applicable |
| **Output / Evidence** | Steps to see the issue; screenshot or log references | When applicable |
| **Priority** | Numeric (e.g., 8/10) with rationale | Yes |
| **Impact** | Bottleneck to x, y, z; who/what is blocked | When applicable |
| **Benefits** | What doing this unlocks (x, y, z) | When applicable |
| **Expected Behavior** | How it should work | Yes |
| **Desired Structure** | The target state we want | When applicable |
| **Acceptance Criteria** | Concrete checks for "done" | Yes |
| **Human Input** | Output of interviewing the human partner; desired outcome | When BLOCKED or human input needed |
| **Questions for Human** | Prompts to ask when info is missing (ARL human-probing) | When BLOCKED |
| **Resources** | Compact table: skills, agents, prior work | Yes |
| **Dependencies** | Items this depends on; items this unblocks | When applicable |
| **Blockers** | Missing prerequisites; RT-ICA BLOCKED reason | When BLOCKED |
| **Effort** | Small / Medium / High | When estimable |

---

## RT-ICA Integration

- **APPROVED**: Item is ready; groomed sections filled
- **BLOCKED**: Include **Blockers** and **Questions for Human**; do not treat as ready for planning

---

## Example Groomed Body

```markdown
## Description

Original problem statement from backlog.

## Groomed (2026-02-23)

### Reproducibility

1. Run `uv run prek run --files foo.py`
2. Observe LK001 on line 42

### Output / Evidence

- Screenshot: `.claude/backlog/assets/p1-foo-screenshot.png`
- Log excerpt in `plugins/plugin-creator/planning/qa-report.md`

### Priority

8/10 — Blocks plugin validation CI gate; affects all contributors.

### Impact

- Blocks: manifest-sync job, pre-commit hooks
- Bottleneck: Every PR runs plugin_validator

### Benefits

- CI gate becomes reliable
- Pre-commit output actionable

### Expected Behavior

Plugin validator should report unique files, not validator invocations.

### Acceptance Criteria

1. Running validator on 1 file shows "Total files: 1"
2. Summary counts unique files

### Resources

| Type | Item |
|------|------|
| Skill | plugin-creator:claude-hooks-reference-2026 |
| Agent | @refactor-executor |
| Prior work | plan/tasks-2-validator-ux-coverage.md |
```
