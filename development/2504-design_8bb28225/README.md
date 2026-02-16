# Design

| Property | Value |
|----------|-------|
| **Name** | Design |
| **Repository** | [rocklambros/zerg](https://raw.githubusercontent.com/rocklambros/zerg/main/.zerg/specs/core-refactoring/design.md) (⭐ 17) |
| **Original Path** | `.zerg/specs/core-refactoring/design.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-06 |
| **Updated** | 2026-02-06 |
| **File Hash** | `8bb28225f8d270f3...` |

## Description

After PR 3 (Registry):

zerg/
├── worker_registry.py        ~150L   WorkerRegistry: threadsafe worker state
└── (5 consumers migrated from raw dict to registry)

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [rocklambros/zerg](https://raw.githubusercontent.com/rocklambros/zerg/main/.zerg/specs/core-refactoring/design.md)*
