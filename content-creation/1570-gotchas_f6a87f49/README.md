# Gotchas

| Property | Value |
|----------|-------|
| **Name** | Gotchas |
| **Repository** | [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/cache-reserve/gotchas.md) (⭐ 392) |
| **Original Path** | `packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/cache-reserve/gotchas.md` |
| **Category** | content-creation |
| **Subcategory** | media |
| **Tags** | content creation |
| **Created** | 2026-02-05 |
| **Updated** | 2026-02-05 |
| **File Hash** | `f6a87f49c1d218fb...` |

## Description

Cause: Asset is not cacheable, TTL < 10 hours, ContentLength header missing, or blocking headers present (SetCookie, Vary: _)  
Solution: Ensure minimum TTL of 10+ hours (CacheControl: public, maxage=36000), add ContentLength header, remove SetCookie header, and set Vary: AcceptEncoding (not _)

**Tags:** `content creation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [tech-leads-club/agent-skills](https://raw.githubusercontent.com/tech-leads-club/agent-skills/main/packages/skills-catalog/skills/(cloud)/cloudflare-deploy/references/cache-reserve/gotchas.md)*
