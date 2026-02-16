# Design

| Property | Value |
|----------|-------|
| **Name** | Design |
| **Repository** | [rocklambros/zerg](https://raw.githubusercontent.com/rocklambros/zerg/main/.gsd/specs/github-issues-batch/design.md) (⭐ 17) |
| **Original Path** | `.gsd/specs/github-issues-batch/design.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-04 |
| **Updated** | 2026-02-04 |
| **File Hash** | `bc198ef7bcf59ab5...` |

## Description

1. Pause signal: Orchestrator sets paused=true in state → Workers check is_paused() before claiming
2. Level transitions: Orchestrator sets current_level → Workers filter tasks by level
3. State validation: All set_task_status() calls validate transition legality

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [rocklambros/zerg](https://raw.githubusercontent.com/rocklambros/zerg/main/.gsd/specs/github-issues-batch/design.md)*
