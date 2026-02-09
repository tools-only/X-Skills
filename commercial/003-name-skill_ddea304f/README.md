# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/databricks-pack/skills/databricks-upgrade-migration/SKILL.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/databricks-pack/skills/databricks-upgrade-migration/SKILL.md` |
| **Category** | commercial |
| **Subcategory** | ecommerce |
| **Tags** | commercial |
| **Created** | 2026-01-26 |
| **Updated** | 2026-01-26 |
| **File Hash** | `ddea304fb0d9bd6e...` |

## Description

def upgrade_cluster_dbr(
    w: WorkspaceClient,
    cluster_id: str,
    target_version: str = "14.3.xscala2.12",
    dry_run: bool = True,
) > dict:
    """
    Upgrade cluster to new DBR version.

**Tags:** `commercial`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/databricks-pack/skills/databricks-upgrade-migration/SKILL.md)*
