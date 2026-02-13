# Tracing

| Property | Value |
|----------|-------|
| **Name** | Tracing |
| **Repository** | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/tracing.md) (⭐ 612) |
| **Original Path** | `kernel/subsystem/tracing.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-11 |
| **Updated** | 2026-02-12 |
| **File Hash** | `516de2962436b656...` |

## Description

Incorrect use of TRACE_EVENT macros causes data corruption in the ring
buffer, truncated strings, or kernel panics from buffer overflows. Misusing
TP_fast_assign() with side effects leads to behavioral differences
depending on whether tracing is enabled, creating hardtodiagnose bugs.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/tracing.md)*
