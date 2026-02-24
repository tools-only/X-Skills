# Task File Format

YAML frontmatter-based task file format for SAM workflows.

---

## Specification

Full specification: [TASK_FILE_FORMAT.md](../../.claude/docs/TASK_FILE_FORMAT.md)

---

## Key Fields

| Field | Type | Required |
|-------|------|----------|
| task | string | Yes |
| title | string | Yes |
| status | enum | Yes |
| agent | string | No |
| dependencies | array | No |
| priority | integer | No |
| complexity | enum | No |

---

## File Organization

- **Single file**: Multi-task with `---` delimiters
- **Directory**: One task per file, `{task-id}-{slug}.md`

---

## Sections

Context, Objective, Requirements, Constraints, Expected Outputs, Acceptance Criteria, Verification Steps, Handoff.
