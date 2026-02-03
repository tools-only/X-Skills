# 2026 02 02 Remove Heuristic Token Counting

| Property | Value |
|----------|-------|
| **Name** | 2026 02 02 Remove Heuristic Token Counting |
| **Repository** | [alchemiststudiosDOTai/tunacode](https://raw.githubusercontent.com/alchemiststudiosDOTai/tunacode/master/.claude/delta/2026-02-02-remove-heuristic-token-counting.md) (⭐ 112) |
| **Original Path** | `.claude/delta/2026-02-02-remove-heuristic-token-counting.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-03 |
| **Updated** | 2026-02-03 |
| **File Hash** | `9c1ea9ff5df2e789...` |

## Description

The agent run flow no longer recomputes token totals by iterating over message
history. This removes the slow heuristic update_token_count path in favor of
usage totals as the canonical source for token display.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [alchemiststudiosDOTai/tunacode](https://raw.githubusercontent.com/alchemiststudiosDOTai/tunacode/master/.claude/delta/2026-02-02-remove-heuristic-token-counting.md)*
