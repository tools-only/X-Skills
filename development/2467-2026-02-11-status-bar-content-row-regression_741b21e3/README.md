# 2026 02 11 Status Bar Content Row Regression

| Property | Value |
|----------|-------|
| **Name** | 2026 02 11 Status Bar Content Row Regression |
| **Repository** | [alchemiststudiosDOTai/tunacode](https://raw.githubusercontent.com/alchemiststudiosDOTai/tunacode/master/.claude/delta/2026-02-11-status-bar-content-row-regression.md) (⭐ 113) |
| **Original Path** | `.claude/delta/2026-02-11-status-bar-content-row-regression.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-11 |
| **Updated** | 2026-02-11 |
| **File Hash** | `741b21e3d7e31795...` |

## Description

The status bar gained a top bevel border during the UI revamp, but its height remained 1. In Textual, border rows consume layout space, so a onerow widget with a top border leaves zero content rows. This caused bottombar information text to render incorrectly / disappear.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [alchemiststudiosDOTai/tunacode](https://raw.githubusercontent.com/alchemiststudiosDOTai/tunacode/master/.claude/delta/2026-02-11-status-bar-content-row-regression.md)*
