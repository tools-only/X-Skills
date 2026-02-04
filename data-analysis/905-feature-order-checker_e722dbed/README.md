# Feature Order Checker

| Property | Value |
|----------|-------|
| **Name** | Feature Order Checker |
| **Repository** | [koreyba/Claude-Skill-pdf-to-epub](https://raw.githubusercontent.com/koreyba/Claude-Skill-pdf-to-epub/master/docs/ai/design/feature-order-checker.md) (⭐ 15) |
| **Original Path** | `docs/ai/design/feature-order-checker.md` |
| **Category** | data-analysis |
| **Subcategory** | visualization |
| **Tags** | data analysis |
| **Created** | 2025-12-30 |
| **Updated** | 2025-12-30 |
| **File Hash** | `e722dbedd7e35527...` |

## Description

mermaid
graph TD
    A[Source PDF Text] >|Segment| B[CompletenessChecker]
    C[Target EPUB Text] >|Scan| B
    B >|Found Chunks + Positions| D[OrderChecker]
    D >|Calculate LIS| E[Order Score]

**Tags:** `data analysis`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [koreyba/Claude-Skill-pdf-to-epub](https://raw.githubusercontent.com/koreyba/Claude-Skill-pdf-to-epub/master/docs/ai/design/feature-order-checker.md)*
