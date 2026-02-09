# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/databricks-pack/skills/databricks-sdk-patterns/SKILL.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/databricks-pack/skills/databricks-sdk-patterns/SKILL.md` |
| **Category** | development |
| **Subcategory** | tools |
| **Tags** | development |
| **Created** | 2026-01-26 |
| **Updated** | 2026-01-26 |
| **File Hash** | `33b2324a9f5415e6...` |

## Description

@lru_cache(maxsize=1)
def get_workspace_client(profile: Optional[str] = None) > WorkspaceClient:
    """Get singleton WorkspaceClient instance."""
    return WorkspaceClient(profile=profile)

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/databricks-pack/skills/databricks-sdk-patterns/SKILL.md)*
