# Js Batch Dom Css

| Property | Value |
|----------|-------|
| **Name** | Js Batch Dom Css |
| **Repository** | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/vercel-react-best-practices/rules/js-batch-dom-css.md) (⭐ 13) |
| **Original Path** | `skills/vercel-react-best-practices/rules/js-batch-dom-css.md` |
| **Category** | automation |
| **Subcategory** | scripting |
| **Tags** | automation |
| **Created** | 2026-01-15 |
| **Updated** | 2026-01-20 |
| **File Hash** | `480b891b9eaf96e9...` |

## Description

Avoid interleaving style writes with layout reads. When you read a layout property (like offsetWidth, getBoundingClientRect(), or getComputedStyle()) between style changes, the browser is forced to trigger a synchronous reflow.

**Tags:** `automation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/vercel-react-best-practices/rules/js-batch-dom-css.md)*
