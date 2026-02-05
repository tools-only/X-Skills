# Patterns

| Property | Value |
|----------|-------|
| **Name** | Patterns |
| **Repository** | [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/bot-management/patterns.md) (⭐ 392) |
| **Original Path** | `packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/bot-management/patterns.md` |
| **Category** | research |
| **Subcategory** | data-gathering |
| **Tags** | research |
| **Created** | 2026-02-05 |
| **Updated** | 2026-02-05 |
| **File Hash** | `d3383b422f68feba...` |

## Description

txt
 High security for checkout
(cf.bot_management.score lt 50 and http.request.uri.path in {"/checkout" "/cart/add"} and not cf.bot_management.verified_bot and not cf.bot_management.corporate_proxy)
Action: Managed Challenge

**Tags:** `research`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/bot-management/patterns.md)*
