# Requirements

| Property | Value |
|----------|-------|
| **Name** | Requirements |
| **Repository** | [rocklambros/zerg](https://raw.githubusercontent.com/rocklambros/zerg/main/.gsd/specs/rush-perf-fix/requirements.md) (⭐ 17) |
| **Original Path** | `.gsd/specs/rush-perf-fix/requirements.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-04 |
| **Updated** | 2026-02-04 |
| **File Hash** | `35e93558ec7f2b89...` |

## Description

ZERG rush execution is PAINFULLY slow due to excessive Docker subprocess calls in the container monitoring loop. The poll loop runs every 5 seconds and executes docker inspect + docker exec for EVERY worker, causing 120+ Docker calls per minute.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [rocklambros/zerg](https://raw.githubusercontent.com/rocklambros/zerg/main/.gsd/specs/rush-perf-fix/requirements.md)*
