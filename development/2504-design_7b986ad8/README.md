# Design

| Property | Value |
|----------|-------|
| **Name** | Design |
| **Repository** | [rocklambros/zerg](https://raw.githubusercontent.com/rocklambros/zerg/main/.gsd/specs/github-issue-98/design.md) (⭐ 17) |
| **Original Path** | `.gsd/specs/github-issue-98/design.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-04 |
| **Updated** | 2026-02-04 |
| **File Hash** | `7b986ad8efe01796...` |

## Description

1. User runs /zerg:design
2. Claude generates taskgraph.json following design.core.md instructions
3. NEW: Instructions mandate appending wiring verification task to L5
4. Wiring task depends on all L4 tasks
5. Workers execute L1L4 tasks
6. Final worker executes wiring task with scoped pytest

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [rocklambros/zerg](https://raw.githubusercontent.com/rocklambros/zerg/main/.gsd/specs/github-issue-98/design.md)*
