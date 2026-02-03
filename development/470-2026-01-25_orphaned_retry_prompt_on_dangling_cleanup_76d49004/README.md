# 2026 01 25 Orphaned Retry Prompt On Dangling Cleanup

| Property | Value |
|----------|-------|
| **Name** | 2026 01 25 Orphaned Retry Prompt On Dangling Cleanup |
| **Repository** | [alchemiststudiosDOTai/tunacode](https://raw.githubusercontent.com/alchemiststudiosDOTai/tunacode/master/.claude/delta/2026-01-25_orphaned_retry_prompt_on_dangling_cleanup.md) (⭐ 112) |
| **Original Path** | `.claude/delta/2026-01-25_orphaned_retry_prompt_on_dangling_cleanup.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-01-25 |
| **Updated** | 2026-01-25 |
| **File Hash** | `76d49004507d3d4b...` |

## Description

Fixed a bug where retryprompt parts were left orphaned in message history after their corresponding toolcall parts were pruned during dangling tool call cleanup. The fix generalizes the filter to remove ANY part with a tool_call_id matching a dangling ID, not just toolcall parts.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [alchemiststudiosDOTai/tunacode](https://raw.githubusercontent.com/alchemiststudiosDOTai/tunacode/master/.claude/delta/2026-01-25_orphaned_retry_prompt_on_dangling_cleanup.md)*
