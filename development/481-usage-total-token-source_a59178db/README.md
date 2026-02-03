# Usage Total Token Source

| Property | Value |
|----------|-------|
| **Name** | Usage Total Token Source |
| **Repository** | [alchemiststudiosDOTai/tunacode](https://raw.githubusercontent.com/alchemiststudiosDOTai/tunacode/master/.claude/patterns/usage-total-token-source.md) (⭐ 112) |
| **Original Path** | `.claude/patterns/usage-total-token-source.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-03 |
| **Updated** | 2026-02-03 |
| **File Hash** | `a59178db08226712...` |

## Description

The canonical token total for UI display is the accumulated pydanticai usage totals
(prompt + completion) from usage.session_total_usage, not heuristic message counting.
This aligns the UI with billed API usage and avoids slow permessage estimation.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [alchemiststudiosDOTai/tunacode](https://raw.githubusercontent.com/alchemiststudiosDOTai/tunacode/master/.claude/patterns/usage-total-token-source.md)*
