# Model Migration

| Property | Value |
|----------|-------|
| **Name** | Model Migration |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/openrouter-pack/skills/openrouter-upgrade-migration/references/model-migration.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/openrouter-pack/skills/openrouter-upgrade-migration/references/model-migration.md` |
| **Category** | content-creation |
| **Subcategory** | writing |
| **Tags** | content creation |
| **Created** | 2026-01-06 |
| **Updated** | 2026-01-06 |
| **File Hash** | `367d3370eef2e2c1...` |

## Description

def migrate_model_name(model: str) > str:
    """Migrate deprecated model to current equivalent."""
    return MODEL_MIGRATIONS.get(model, model)

**Tags:** `content creation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/openrouter-pack/skills/openrouter-upgrade-migration/references/model-migration.md)*
