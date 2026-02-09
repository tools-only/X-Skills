# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/langchain-pack/skills/langchain-rate-limits/SKILL.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/langchain-pack/skills/langchain-rate-limits/SKILL.md` |
| **Category** | automation |
| **Subcategory** | scripting |
| **Tags** | automation |
| **Created** | 2026-01-09 |
| **Updated** | 2026-01-09 |
| **File Hash** | `b9e76aafe7d4cf14...` |

## Description

@retry(
    stop=stop_after_attempt(5),
    wait=wait_exponential(multiplier=1, min=4, max=60),
    retry=retry_if_exception_type((RateLimitError, APIError))
)
def call_with_retry(chain, input_data):
    """Call chain with exponential backoff."""
    return chain.invoke(input_data)

**Tags:** `automation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/langchain-pack/skills/langchain-rate-limits/SKILL.md)*
