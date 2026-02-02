# Equation Counting Bug Fix

| Property | Value |
|----------|-------|
| **Name** | Equation Counting Bug Fix |
| **Repository** | [dmccreary/claude-skills](https://raw.githubusercontent.com/dmccreary/claude-skills/main/src/book-metrics/EQUATION_COUNT_FIX.md) (⭐ 14) |
| **Original Path** | `src/book-metrics/EQUATION_COUNT_FIX.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2025-11-25 |
| **Updated** | 2025-11-25 |
| **File Hash** | `e804b0bc980a21c0...` |

## Description

The count_equations_in_file() function in bookmetrics.py was significantly overcounting LaTeX equations. For example, a signal processing textbook showed 2,421 equations when the actual count was much lower.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [dmccreary/claude-skills](https://raw.githubusercontent.com/dmccreary/claude-skills/main/src/book-metrics/EQUATION_COUNT_FIX.md)*
