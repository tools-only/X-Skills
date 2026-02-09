# Gotchas

| Property | Value |
|----------|-------|
| **Name** | Gotchas |
| **Repository** | [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/cache-reserve/gotchas.md) (🔥 19.8k) |
| **Original Path** | `cli-tool/components/skills/development/cloudflare-deploy/references/cache-reserve/gotchas.md` |
| **Category** | content-creation |
| **Subcategory** | media |
| **Tags** | content creation |
| **Created** | 2026-02-08 |
| **Updated** | 2026-02-08 |
| **File Hash** | `3642800cd57f39fa...` |

## Description

Cause: Asset is not cacheable, TTL < 10 hours, ContentLength header missing, or blocking headers present (SetCookie, Vary: )  
Solution: Ensure minimum TTL of 10+ hours (CacheControl: public, maxage=36000), add ContentLength header, remove SetCookie header, and set Vary: AcceptEncoding (not )

**Tags:** `content creation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/cli-tool/components/skills/development/cloudflare-deploy/references/cache-reserve/gotchas.md)*
