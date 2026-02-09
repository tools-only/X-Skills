# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/databricks-pack/skills/databricks-rate-limits/SKILL.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/databricks-pack/skills/databricks-rate-limits/SKILL.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-01-26 |
| **Updated** | 2026-01-26 |
| **File Hash** | `ad4e2c1177653127...` |

## Description

def with_exponential_backoff(
    operation: Callable[[], T],
    max_retries: int = 5,
    base_delay: float = 1.0,
    max_delay: float = 60.0,
    jitter_factor: float = 0.5,
) > T:
    """
    Execute operation with exponential backoff on rate limits.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/databricks-pack/skills/databricks-rate-limits/SKILL.md)*
