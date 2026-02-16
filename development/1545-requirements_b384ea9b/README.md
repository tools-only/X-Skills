# Requirements

| Property | Value |
|----------|-------|
| **Name** | Requirements |
| **Repository** | [rocklambros/zerg](https://raw.githubusercontent.com/rocklambros/zerg/main/.gsd/specs/resolving-the-rest-of-134-completely/requirements.md) (⭐ 17) |
| **Original Path** | `.gsd/specs/resolving-the-rest-of-134-completely/requirements.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-15 |
| **Updated** | 2026-02-15 |
| **File Hash** | `b384ea9bc1aa1afb...` |

## Description

Issue 134 identified that ZERG performs redundant rglob directory traversals across multiple modules. A shared fs_utils.collect_files() utility was created (PR 161) and adopted by repo_map.py and stack_detector.py, but 18 rglob calls across 16 files remain unconverted.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [rocklambros/zerg](https://raw.githubusercontent.com/rocklambros/zerg/main/.gsd/specs/resolving-the-rest-of-134-completely/requirements.md)*
