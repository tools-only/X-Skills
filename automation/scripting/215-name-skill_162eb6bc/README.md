# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/databricks-pack/skills/databricks-performance-tuning/SKILL.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/databricks-pack/skills/databricks-performance-tuning/SKILL.md` |
| **Category** | automation |
| **Subcategory** | scripting |
| **Tags** | automation |
| **Created** | 2026-01-26 |
| **Updated** | 2026-01-26 |
| **File Hash** | `162eb6bc1eb61116...` |

## Description

python
 Cluster sizing calculator
def recommend_cluster_size(
    data_size_gb: float,
    complexity: str = "medium",   low, medium, high
    parallelism_need: str = "standard",   standard, high
) > dict:
    """
    Recommend cluster configuration based on workload.

**Tags:** `automation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/databricks-pack/skills/databricks-performance-tuning/SKILL.md)*
