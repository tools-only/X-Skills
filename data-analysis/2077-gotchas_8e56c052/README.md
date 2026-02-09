# Gotchas

| Property | Value |
|----------|-------|
| **Name** | Gotchas |
| **Repository** | [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/smart-placement/gotchas.md) (🔥 19.8k) |
| **Original Path** | `cli-tool/components/skills/development/cloudflare-deploy/references/smart-placement/gotchas.md` |
| **Category** | data-analysis |
| **Subcategory** | processing |
| **Tags** | data analysis |
| **Created** | 2026-02-08 |
| **Updated** | 2026-02-08 |
| **File Hash** | `8e56c052e6b513b9...` |

## Description

Cause: Not enough traffic for Smart Placement to analyze
Solution:
 Ensure Worker receives consistent global traffic
 Wait longer (analysis takes up to 15 minutes)
 Send test traffic from multiple global locations
 Check Worker has fetch event handler

**Tags:** `data analysis`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/smart-placement/gotchas.md)*
