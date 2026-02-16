# Design

| Property | Value |
|----------|-------|
| **Name** | Design |
| **Repository** | [rocklambros/zerg](https://raw.githubusercontent.com/rocklambros/zerg/main/.gsd/specs/fix-task-list-id/design.md) (⭐ 17) |
| **Original Path** | `.gsd/specs/fix-task-list-id/design.md` |
| **Category** | daily-assistant |
| **Subcategory** | tasks |
| **Tags** | daily assistant |
| **Created** | 2026-01-31 |
| **Updated** | 2026-01-31 |
| **File Hash** | `b2066ea21da90433...` |

## Description

Orchestrator (no env var set)  →  uses DEFAULT task list
    │
    └── spawns worker with CLAUDE_CODE_TASK_LIST_ID=feature
            → uses FEATURE task list  ← MISMATCH

**Tags:** `daily assistant`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [rocklambros/zerg](https://raw.githubusercontent.com/rocklambros/zerg/main/.gsd/specs/fix-task-list-id/design.md)*
