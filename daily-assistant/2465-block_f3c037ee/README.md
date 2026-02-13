# Block

| Property | Value |
|----------|-------|
| **Name** | Block |
| **Repository** | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/block.md) (⭐ 612) |
| **Original Path** | `kernel/subsystem/block.md` |
| **Category** | daily-assistant |
| **Subcategory** | notes |
| **Tags** | daily assistant |
| **Created** | 2026-02-11 |
| **Updated** | 2026-02-12 |
| **File Hash** | `f3c037ee8a631a27...` |

## Description

Failing to hold a queue freeze during teardown or reconfiguration allows bios
to complete concurrently, causing useafterfree on queue state, stale
elevator data, or torn reads of q>nr_requests.

**Tags:** `daily assistant`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/block.md)*
