# Reduce Sdk Overhead

| Property | Value |
|----------|-------|
| **Name** | Reduce Sdk Overhead |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/sentry-pack/skills/sentry-performance-tuning/references/reduce-sdk-overhead.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/sentry-pack/skills/sentry-performance-tuning/references/reduce-sdk-overhead.md` |
| **Category** | automation |
| **Subcategory** | scripting |
| **Tags** | automation |
| **Created** | 2026-01-06 |
| **Updated** | 2026-01-06 |
| **File Hash** | `2edc86a60c765e41...` |

## Description

// Only include needed integrations
  integrations: [
    new Sentry.Integrations.Http({ tracing: true }),
    // Remove unused integrations to reduce overhead
  ],

**Tags:** `automation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/sentry-pack/skills/sentry-performance-tuning/references/reduce-sdk-overhead.md)*
